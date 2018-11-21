/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    //TransferFunction2D tFunc2D;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    
    // Function types
    public static final int SLICER = 1;
    public static final int MIP = 2;
    public static final int COMPOSITING = 3;
    public static final int FUNC2D = 4;
    public static final int MAXGRADIENT = 5;

    private int function;
    
    final static int INTERACTIVE_MODE_STEP = 20;
    final static int INTERACTIVE_MODE_GRANULARITY = 2;
    final static int NON_INTERACTIVE_MODE_STEP = 1;
    final static int NON_INTERACTIVE_MODE_GRANULARITY = 1;

    // Ray parameters
    int granularity = NON_INTERACTIVE_MODE_GRANULARITY;
    int step = NON_INTERACTIVE_MODE_STEP;
    
    // Shading option
    private boolean shading;
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        tFunc.setTestFunc();
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }
    
    public void setShading(boolean shading) {
        this.shading = shading;
    }
    
    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
    
    // This function does a nearest neighbor interpolation
    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return volume.getVoxel(x, y, z);
    }
    
    private void clearImage(){
    	for (int j = 0; j < image.getHeight(); j++) {
	        for (int i = 0; i < image.getWidth(); i++) {
	            image.setRGB(i, j, 0);
	        }
	    }
    }

    void slicer(double[] viewMatrix) {

        // clear image
        clearImage();

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        // Iterate on every pixel
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                int val = volume.getVoxelLinearInterpolate(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
                
                setRGB2Image(voxelColor, i, j, granularity);
            }
        }

    }
    
    void setRGB2Image(TFColor voxelColor, int i, int j, int granularity) { //
	    // BufferedImage expects a pixel color packed as ARGB in an int 
	    int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
	    int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
	    int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
	    int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
	    int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
	    for (int m = 0; m < granularity; m++) {
	        for (int n = 0; n < granularity; n++) {
	              if (n + i >= image.getWidth() || n + i < 0 || m + j >= image.getWidth() || m + j < 0) {
	              } else {
	                  image.setRGB(n + i, m + j, pixelColor);
	            }
	        }
	    }
    }
    
    void raycast(double[] viewMatrix) {

    	 // clear image
        clearImage();
        
        if(this.interactiveMode){
            granularity = INTERACTIVE_MODE_GRANULARITY; // trying to make it faster during interaction
        }else {
            granularity = NON_INTERACTIVE_MODE_GRANULARITY;
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] entryPoint = new double[3];
        double[] exitPoint = new double[3];
        
        //The ray is pointing towards the scene
        double[] rayVector = new double[3];
        rayVector[0]=-viewVec[0];
        rayVector[1]=-viewVec[1];
        rayVector[2]=-viewVec[2];
        
        // We use orthographic projection. Viewer is far away at the infinite, all pixels have the same rayVector.
        
        // ray computation for each pixel
        for (int j = 0; j < image.getHeight(); j += granularity) {
            for (int i = 0; i < image.getWidth(); i += granularity) {
                // compute starting points of rays in a plane shifted backwards to a position behind the data set
            	
            	// Pixel coordinate is calculate having the center (0,0) of the view plane aligned with the center of the volume and moved a distance equivalent
                // to the diagonal to make sure I am far away enough.
                double diagonal = Math.sqrt((volume.getDimX()*volume.getDimX())+(volume.getDimY()*volume.getDimY())+ (volume.getDimZ()*volume.getDimZ()))/2;               
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * diagonal + volume.getDimX() / 2.0;
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * diagonal + volume.getDimY() / 2.0;
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * diagonal + volume.getDimZ() / 2.0;
            	
            	
            	// compute the entry and exit point of the ray
                computeEntryAndExit(pixelCoord, rayVector, entryPoint, exitPoint);
                if ((entryPoint[0] > -1.0) && (exitPoint[0] > -1.0)) {
                    int val = 0;
                    val = traceRayMIP(entryPoint, exitPoint, rayVector, step);
                    for (int ii = i; ii < i + granularity; ii++) {
                        for (int jj = j; jj < j + granularity; jj++) {
                            image.setRGB(ii, jj, val);
                        }
                    }
                }

            }
        }
    }
    
    void composite(double[] viewMatrix) {

   	   // clear image
       clearImage();
       
       if(this.interactiveMode){
           granularity = INTERACTIVE_MODE_GRANULARITY; // trying to make it faster during interaction
       }else {
           granularity = NON_INTERACTIVE_MODE_GRANULARITY;
       }

       // vector uVec and vVec define a plane through the origin, 
       // perpendicular to the view vector viewVec
       double[] viewVec = new double[3];
       double[] uVec = new double[3];
       double[] vVec = new double[3];
       VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
       VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
       VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

       // image is square
       int imageCenter = image.getWidth() / 2;

       double[] pixelCoord = new double[3];
       double[] entryPoint = new double[3];
       double[] exitPoint = new double[3];
       
       //The ray is pointing towards the scene
       double[] rayVector = new double[3];
       rayVector[0]=-viewVec[0];
       rayVector[1]=-viewVec[1];
       rayVector[2]=-viewVec[2];
       
       // We use orthographic projection. Viewer is far away at the infinite, all pixels have the same rayVector.
       
       // ray computation for each pixel
       for (int j = 0; j < image.getHeight(); j += granularity) {
           for (int i = 0; i < image.getWidth(); i += granularity) {
               // compute starting points of rays in a plane shifted backwards to a position behind the data set
           	computePixelCoordinatesBehind(pixelCoord,viewVec,uVec,vVec,i,j);
           	// compute the entry and exit point of the ray
               computeEntryAndExit(pixelCoord, rayVector, entryPoint, exitPoint);
               if ((entryPoint[0] > -1.0) && (exitPoint[0] > -1.0)) {
                   int val = 0;
                   val = traceRayComposite(entryPoint, exitPoint, rayVector, step);
                   for (int ii = i; ii < i + granularity; ii++) {
                       for (int jj = j; jj < j + granularity; jj++) {
                           image.setRGB(ii, jj, val);
                       }
                   }
               }

           }
       }
    }
    
    //Implementation of the MIP per ray  given the entry and exit point and the ray direction
    // sampleStep indicates the distance between samples
    int traceRayMIP(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep) {
    	//compute the increment and the number of samples you need and iterate over them.
 
    	// fast mode: speed up rendering by decrease resolution
        if(this.interactiveMode) {
            sampleStep *= INTERACTIVE_MODE_STEP;
        } else {
            sampleStep *= NON_INTERACTIVE_MODE_STEP;
        }
        
        //get volume max 
        double max=volume.getMaximum();
        
        //initializing the iteration from the entry point
        double pos[] = new double[3];
        pos[0] = entryPoint[0];
        pos[1] = entryPoint[1];
        pos[2] = entryPoint[2];
        
        //calculating number of steps
        double totalDistance = VectorMath.distance(entryPoint, exitPoint);
        int numStep = (int)Math.floor(totalDistance/sampleStep);
        // initializing intensities for this ray
        double maxIntensity = 0;
        double currIntensity = 0;
        
        //iterating front to back
        for(int steps = 0; steps<numStep; steps++){
            
        	// speed up rendering when moving volume 

            if(this.interactiveMode){
                currIntensity = getVoxel(pos); // nearest neighbor calculation during interaction
            }else {
                currIntensity = volume.getVoxelLinearInterpolate(pos);
            }
            
            //store maximum value
            if(currIntensity>maxIntensity){
                maxIntensity = currIntensity;   
            }
            
            //update position
            pos[0] += rayVector[0]*sampleStep;
            pos[1] += rayVector[1]*sampleStep;
            pos[2] += rayVector[2]*sampleStep;
        }
        // Map the intensity to a grey value by linear scaling
        double r = maxIntensity/max;
        double g = r;
        double b = r;
        //if the maximum = 0 make the voxel transparent
        double alpha = maxIntensity > 0 ? 1.0: 0.0 ;
                
        int color = computeImageColor(r,g,b,alpha);
        return color;
    }
    
    int traceRayComposite(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep) {

    	// fast mode: speed up rendering by decrease resolution
        if(this.interactiveMode) {
            sampleStep *= INTERACTIVE_MODE_STEP;
        } else {
            sampleStep *= NON_INTERACTIVE_MODE_STEP;
        }

        //initializing the iteration from the entry point
        double pos[] = new double[3];
        pos[0] = entryPoint[0];
        pos[1] = entryPoint[1];
        pos[2] = entryPoint[2];
        
        //calculating number of steps
        double totalDistance = VectorMath.distance(entryPoint, exitPoint);
        int numStep = (int)Math.floor(totalDistance/sampleStep);
        
        // initializing colors for this pixel
        double currIntensity = 0;
        TFColor finalColor = new TFColor(0,0,0,0);
        TFColor tf2dColor = new TFColor(0,0,0,0);
        double r = 0;
        double g = 0;
        double b = 0;
        double a = 0;
          
        VoxelGradient gradient = new VoxelGradient();
        
        //iterating front to back
        for(int steps = 0; steps<numStep; steps++){
            // 1D transfer function
         
            // speed up rendering when moving volume 

            if(this.interactiveMode){
                currIntensity = getVoxel(pos); // nearest neighbor calculation during interaction
            }else {
                currIntensity = volume.getVoxelLinearInterpolate(pos);
            }
            TFColor  currColor = tFunc.getColor((int)currIntensity);
            
            //updating color for front to back compositing
            r = currColor.r;
            g = currColor.g;
            b = currColor.b;
            a = currColor.a;
             
            //update position
            pos[0] += rayVector[0]*sampleStep;
            pos[1] += rayVector[1]*sampleStep;
            pos[2] += rayVector[2]*sampleStep;
            

            //color compositing - front to back
            finalColor.r = finalColor.r + (1 - finalColor.a)*a*r;
            finalColor.g = finalColor.g + (1 - finalColor.a)*a*g;
            finalColor.b = finalColor.b + (1 - finalColor.a)*a*b;
            finalColor.a = finalColor.a + (1 - finalColor.a)*a;
            
            //early ray termination
            if(finalColor.a>=0.98){
                break;
            }

        }

        double rn = finalColor.r;
        double gn = finalColor.g;
        double bn = finalColor.b;
        double alpha = finalColor.a ;  
        
        int color = computeImageColor(rn,gn,bn,alpha);
        return color;
    }
    
    int computeImageColor(double r, double g, double b, double a){
		int c_alpha = 	a <= 1.0 ? (int) Math.floor(a * 255) : 255;
        int c_red = 	r <= 1.0 ? (int) Math.floor(r * 255) : 255;
        int c_green = 	g <= 1.0 ? (int) Math.floor(g * 255) : 255; 
        int c_blue = 	b <= 1.0 ? (int) Math.floor(b * 255) : 255;
        int pixelColor = getColorInteger(c_red,c_green,c_blue,c_alpha);
        return pixelColor;
	}
    
    public int getColorInteger(int c_red, int c_green, int c_blue, int c_alpha) {
    	int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
    	return pixelColor;
    } 
    
    void computeEntryAndExit(double[] p, double[] viewVec, double[] entryPoint, double[] exitPoint) {

        for (int i = 0; i < 3; i++) {
            entryPoint[i] = -1;
            exitPoint[i] = -1;
        }

        double[] plane_pos = new double[3];
        double[] plane_normal = new double[3];
        double[] intersection = new double[3];

        VectorMath.setVector(plane_pos, volume.getDimX(), 0, 0);
        VectorMath.setVector(plane_normal, 1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, -1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, volume.getDimY(), 0);
        VectorMath.setVector(plane_normal, 0, 1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, -1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, volume.getDimZ());
        VectorMath.setVector(plane_normal, 0, 0, 1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, 0, -1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

    }
    
    private void intersectFace(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection,
            double[] entryPoint, double[] exitPoint) {

        boolean intersect = intersectLinePlane(plane_pos, plane_normal, line_pos, line_dir,
                intersection);
        if (intersect) {

            double xpos0 = 0;
            double xpos1 = volume.getDimX();
            double ypos0 = 0;
            double ypos1 = volume.getDimY();
            double zpos0 = 0;
            double zpos1 = volume.getDimZ();

            if (validIntersection(intersection, xpos0, xpos1, ypos0, ypos1,
                    zpos0, zpos1)) {
                if (VectorMath.dotproduct(line_dir, plane_normal) < 0) {
                    entryPoint[0] = intersection[0];
                    entryPoint[1] = intersection[1];
                    entryPoint[2] = intersection[2];
                } else {
                    exitPoint[0] = intersection[0];
                    exitPoint[1] = intersection[1];
                    exitPoint[2] = intersection[2];
                }
            }
        }
    }
    
    private boolean validIntersection(double[] intersection, double xb, double xe, double yb,
            double ye, double zb, double ze) {

        return (((xb - 0.5) <= intersection[0]) && (intersection[0] <= (xe + 0.5))
                && ((yb - 0.5) <= intersection[1]) && (intersection[1] <= (ye + 0.5))
                && ((zb - 0.5) <= intersection[2]) && (intersection[2] <= (ze + 0.5)));

    }
    
    private boolean intersectLinePlane(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection) {

        double[] tmp = new double[3];

        for (int i = 0; i < 3; i++) {
            tmp[i] = plane_pos[i] - line_pos[i];
        }

        double denom = VectorMath.dotproduct(line_dir, plane_normal);
        if (Math.abs(denom) < 1.0e-8) {
            return false;
        }

        double t = VectorMath.dotproduct(tmp, plane_normal) / denom;

        for (int i = 0; i < 3; i++) {
            intersection[i] = line_pos[i] + t * line_dir[i];
        }

        return true;
    }

    void computePixelCoordinatesBehind(double pixelCoord[], double viewVec[], double uVec[], double vVec[], int i, int j) {
		int imageCenter = image.getWidth()/2;
                // Pixel coordinate is calculate having the center (0,0) of the view plane aligned with the center of the volume and moved a distance equivalent
                // to the diagonal to make sure I am far away enough.
                double diagonal = Math.sqrt((volume.getDimX()*volume.getDimX())+(volume.getDimY()*volume.getDimY())+ (volume.getDimZ()*volume.getDimZ()))/2;               
		pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * diagonal + volume.getDimX() / 2.0;
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * diagonal + volume.getDimY() / 2.0;
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * diagonal + volume.getDimZ() / 2.0;
	}

    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    public void setFunction(int function) {
        this.function = function;
    }


    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        
        step = this.interactiveMode ? INTERACTIVE_MODE_STEP : NON_INTERACTIVE_MODE_STEP;
        granularity = this.interactiveMode ? INTERACTIVE_MODE_GRANULARITY : NON_INTERACTIVE_MODE_GRANULARITY;
        
        switch (function) {
	        case SLICER:
	            slicer(viewMatrix);
	            break;
	        case MIP:
	        	raycast(viewMatrix);
	            break;
	        case COMPOSITING:
	            
	            break;
	        case FUNC2D:
	            break;
	        case MAXGRADIENT:
	            break;
	        default:
	            slicer(viewMatrix);
	            break;
        }  
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (TFChangeListener listener : listeners) {
            listener.changed();
        }
    }
}
