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
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    
    // Function types
    public static final int SLICER = 1;
    public static final int MIP = 2;
    public static final int COMPOSITING = 3;
    public static final int TRANSFER_FUNCTION_2D = 4;

    private int function;
    
    // Ray step parameter
    final static int INTERACTIVE_MODE_STEP = 4;
    final static int NON_INTERACTIVE_MODE_STEP = 1;
    
    int step = NON_INTERACTIVE_MODE_STEP;
    
    // Resolution parameter
    final static int INTERACTIVE_MODE_RESOLUTION = 2;
    final static int NON_INTERACTIVE_MODE_RESOLUTION = 1;
    
    int resolution = NON_INTERACTIVE_MODE_RESOLUTION;
    
    // Shading option (Phong model implemented)
    private boolean shading = false;
    
    // Phong parameters 
    private final double kAmbient = 0.1;
    private final double kDiffuse = 0.7;
    private final double kSpecular = 0.2;
    private final double kAlpha = 10;
    
    public RaycastRenderer() {
        
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // Set up image for storing the resulting rendering the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY() + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // Create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // Uncomment this to initialize the TF with good starting values for the orange dataset 
        tFunc.setTestFunc();
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
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
    private short getVoxel(double[] coord) {

        // Making sure that point coordinates are within the limits
        if (coord[0] < 0 || coord[0] >= volume.getDimX() || coord[1] < 0 || coord[1] >= volume.getDimY()
                || coord[2] < 0 || coord[2] >= volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return volume.getVoxel(x, y, z);
    }
    
    // This function linearly interpolates the value p0 and p1 given the factor alpha between 0 and 1 
    private float linearInterpolate(float p0, float p1, float alpha) {
        
    	return p0 * (1 - alpha) + p1 * alpha; 
    }
    
    // This function computes the intensity value for the point with coordinates coord, using trilinear interpolation
    private short getVoxelTrilinearInterpolated(double[] coord) {
        
	// Making sure that point coordinates are within the limits.
	if (coord[0] < 0 || coord[0] >= volume.getDimX()-1 || coord[1] < 0 || coord[1] >= volume.getDimY()-1
	        || coord[2] < 0 || coord[2] >= volume.getDimZ()-1) {
	    return 0;
	}
        
	/* We assume that the distance between neighbouring voxels is 1 in all directions. */
	    
	// Identify the cell which contains the point with coordinates coord.
	int x0 = (int) Math.floor(coord[0]); 
	int y0 = (int) Math.floor(coord[1]);
	int z0 = (int) Math.floor(coord[2]);
	    
	/* On a periodic and cubic lattice, let xd, yd and zd be the differences between each of x, y, z and the smaller coordinate 
        related, that is: */
        float xd = (float) (coord[0] - x0);
	float yd = (float) (coord[1] - y0);
	float zd = (float) (coord[2] - z0);
        /* where x0 indicates the lattice point below x, and x1 indicates the lattice point above  x and similarly for 
        y0, y1, z0 and z1. 
        Since we assume the distances between neighbouring voxels are 1 in all directions: 
        the factors (between 0 and 1) needed for interpolate function are (coord-x0)/(x1-x0) = coord-x0 because x1-x0=1. */
	  
	/* Performing trilinear interpolation */
        
        // First we interpolate along x
	float c00 = linearInterpolate(volume.getVoxel(x0, y0, z0), volume.getVoxel(x0 + 1, y0, z0), xd);
	float c10 = linearInterpolate(volume.getVoxel(x0, y0 + 1, z0), volume.getVoxel(x0 + 1, y0 + 1, z0), xd);
	float c01 = linearInterpolate(volume.getVoxel(x0, y0, z0 + 1), volume.getVoxel(x0 + 1, y0, z0 + 1), xd);
	float c11 = linearInterpolate(volume.getVoxel(x0, y0 + 1,z0 + 1), volume.getVoxel(x0 + 1, y0 + 1, z0 + 1), xd); 
	    
	// Interpolate along y 
	float c0 = linearInterpolate(c00, c10, yd);
	float c1 = linearInterpolate(c01, c11, yd);
	    
	// Finally we interpolate these values along z (walking through a line).
	float c = linearInterpolate(c0, c1, zd);
	
        // This gives us a predicted value for the point.   
	return (short) c; 
    }

    // Clear the image
    private void clearImage(){
        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
    }

    /* Slicer function */
    public void slicer(double[] viewMatrix) {

        // Clear image
        this.clearImage();

        // Vector uVec and vVec define a plane through the origin, perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // Image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // Sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];

                // Speed up rendering when moving volume
                int val = this.interactiveMode ? getVoxel(pixelCoord) : getVoxelTrilinearInterpolated(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val / max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                //voxelColor = tFunc.getColor(val);
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                
                image.setRGB(i, j, pixelColor);
            }
        }
    }

    /* Maximum Intensity Projection (MIP) function */
    public void mip(double[] viewMatrix) {

        // Clear image
        this.clearImage();

        // Vector uVec and vVec define a plane through the origin, perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // Image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // Sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        // Step size of samples along the viewing ray
        step = this.interactiveMode ? INTERACTIVE_MODE_STEP : NON_INTERACTIVE_MODE_STEP;
        
        // Ray computation for each pixel
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {

                // Initialize max intensity
                int maxIntensity = 0;
                
                // Steps along the ray to compute the maximum intensity projection for the current pixel
                for (int k = 0; k < volume.getDiagonal(); k+=step) {
             
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * (k - imageCenter) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * (k - imageCenter) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * (k - imageCenter) + volumeCenter[2];

                    // Speed up rendering when moving volume
                    int val = this.interactiveMode ? getVoxel(pixelCoord) : getVoxelTrilinearInterpolated(pixelCoord);
                    
                    // Store maximum value
                    if (val > maxIntensity) {
                        
                        maxIntensity = val;
                    }
                }

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = maxIntensity / max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = maxIntensity > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
        
                image.setRGB(i, j, pixelColor);
            }
        }
    }
    
    /* Compositing ray function */
    void compositing(double[] viewMatrix) {
        
        // Clear the image
        clearImage();

        // Vector uVec and vVec define a plane through the origin, perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        
        // Image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // Sample on a plane through the origin of the volume data
        TFColor voxelColor = new TFColor();
        
        // Step size of samples akong the viewing ray
        step = this.interactiveMode ? INTERACTIVE_MODE_STEP : NON_INTERACTIVE_MODE_STEP;
        
        // Resolution
        resolution = this.interactiveMode ? INTERACTIVE_MODE_RESOLUTION : NON_INTERACTIVE_MODE_RESOLUTION;

        // Ray computation for each pixel
        for (int j = 0; j < image.getHeight(); j += resolution) {
            for (int i = 0; i < image.getWidth(); i += resolution) {
        
                voxelColor = new TFColor(0, 0, 0, 0);
                
                /* Steps along the ray for the current pixel */
                for (int k = 0; k < volume.getDiagonal(); k+=step) {
                    
                    // Get calculate new volumeCenter
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * (k - imageCenter) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * (k - imageCenter) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * (k - imageCenter) + volumeCenter[2];

                    // Speed up rendering when moving volume
                    int val = this.interactiveMode ? getVoxel(pixelCoord) : getVoxelTrilinearInterpolated(pixelCoord);

                    // Get the colors of the Transfer Function
                    TFColor color = tFunc.getColor(val);
                    
                    // Calculate the color of the voxel using the color of the previous voxels (back-to-front compositing order)
                    voxelColor.r = (1 - color.a) * voxelColor.r + color.a*color.r;
                    voxelColor.g = (1 - color.a) * voxelColor.g + color.a*color.g;
                    voxelColor.b = (1 - color.a) * voxelColor.b + color.a*color.b;
                    voxelColor.a = (1 - color.a) * voxelColor.a + color.a;
                }
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                image.setRGB(i, j, pixelColor);
                
                if (this.interactiveMode) {
                    
                    image.setRGB(i, j + 1, pixelColor);
                    image.setRGB(i + 1, j, pixelColor);
                    image.setRGB(i + 1, j + 1, pixelColor);
                }
            }
        }
    }
     
    /* 2D transfer function */
    void transferFunction2D(double[] viewMatrix) {

        // Clear image
        this.clearImage();
        
        // Vector uVec and vVec define a plane through the origin, perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        double[] reverseView = new double[3];
        VectorMath.setVector(reverseView, -viewMatrix[2], -viewMatrix[6], -viewMatrix[10]);

        // Image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        TFColor voxelColor = new TFColor(0, 0, 0, 0);
        double[] rgb;
        
        double alphaV = tfEditor2D.triangleWidget.color.a;
        short baseIntensity = tfEditor2D.triangleWidget.baseIntensity;
        double radius = tfEditor2D.triangleWidget.radius;

        // Thicken the borders in interactive mode
        if(this.interactiveMode){
            
            radius = 2*radius;
        }
        
        // Step size of samples along the viewing ray
        step = this.interactiveMode ? INTERACTIVE_MODE_STEP : NON_INTERACTIVE_MODE_STEP;
        
        // Resolution
        resolution = this.interactiveMode ? INTERACTIVE_MODE_RESOLUTION : NON_INTERACTIVE_MODE_RESOLUTION;
        
        // 2D transfer function computation for each pixel
        for (int j = 0; j < image.getHeight(); j += resolution) {
            for (int i = 0; i < image.getWidth(); i += resolution) {
              
                TFColor compositingColor = new TFColor(0, 0, 0, 0); 
                
                // Steps along the ray for the current pixel
                for (int k = volume.getDiagonal(); k > 0; k-=step) {

                    // Get calculate new volumeCenter
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * (k - imageCenter) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * (k - imageCenter) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * (k - imageCenter) + volumeCenter[2];

                    // Speed up rendering when moving volume
                    int val = this.interactiveMode ? getVoxel(pixelCoord) : getVoxelTrilinearInterpolated(pixelCoord);

                    if (val == 0) {
                        
                        continue;
                    }
                    
                    voxelColor = tfEditor2D.triangleWidget.color.clone();

                    VoxelGradient gradient = gradients.getTriLinearGradient((float) pixelCoord[0], (float) pixelCoord[1], (float) pixelCoord[2]);
                    double opacity = 0;
                    double gradientMag = this.getGradientMag(pixelCoord, gradient);
                    double minGradient = tfEditor2D.triangleWidget.minGradient;
                    double maxGradient = tfEditor2D.triangleWidget.maxGradient;

                    /* Compute the opacity for the 2D transfer function */
                    if (gradientMag <= maxGradient && gradientMag >= minGradient) {
                        
                        if (gradientMag == 0 && val == baseIntensity) {

                            opacity = alphaV;
                        } else if (gradientMag > 0 && (val - radius * gradientMag <= baseIntensity) && (baseIntensity <= val + radius * gradientMag)) {

                            opacity = alphaV * (1 - (1 / radius) * Math.abs((baseIntensity - val) / (gradientMag)));
                        } else {
                            
                            opacity = 0;
                        }
                    } else {
                        
                        opacity = 0;
                    }
                    
                    if (this.shading) { // Shading option (Phong model)
                        
                        double dotProduct = VectorMath.dotproduct(viewVec, gradient.getNormalisedVoxelGradient());
                        rgb = new double[]{0, 0, 0};
                        
                        if (gradientMag > 0 && dotProduct > 0 && opacity > 0) {
                            
                            double[] compRGB = new double[]{voxelColor.r, voxelColor.g, voxelColor.b};
                            
                            for (int z = 0; z < 3; z++) {
                                
                                rgb[z] = this.kAmbient + compRGB[z] * this.kDiffuse * dotProduct + this.kSpecular * Math.pow(dotProduct, this.kAlpha);
                            }
                            
                            voxelColor.r = rgb[0];
                            voxelColor.g = rgb[1];
                            voxelColor.b = rgb[2];
                        } else {
                            
                            continue;
                        }
                    }
                   
                    compositingColor.r = voxelColor.r * opacity + (1 - opacity) * compositingColor.r;
                    compositingColor.g = voxelColor.g * opacity + (1 - opacity) * compositingColor.g;
                    compositingColor.b = voxelColor.b * opacity + (1 - opacity) * compositingColor.b;

                    compositingColor.a = (1 - opacity) * compositingColor.a;
                }

                compositingColor.a = 1 - compositingColor.a;

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = compositingColor.a <= 1.0 ? (int) Math.floor(compositingColor.a * 255) : 255;
                int c_red = compositingColor.r <= 1.0 ? (int) Math.floor(compositingColor.r * 255) : 255;
                int c_green = compositingColor.g <= 1.0 ? (int) Math.floor(compositingColor.g * 255) : 255;
                int c_blue = compositingColor.b <= 1.0 ? (int) Math.floor(compositingColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                image.setRGB(i, j, pixelColor);
                
                if (this.interactiveMode) {
                    
                    image.setRGB(i, j + 1, pixelColor);
                    image.setRGB(i + 1, j, pixelColor);
                    image.setRGB(i + 1, j + 1, pixelColor);
                }
            }
        }
    }
    
    private double getGradientMag(double[] pixelCoord, VoxelGradient gradient) {
        
        if (pixelCoord[0] < 0 || pixelCoord[0] > volume.getDimX() - 1 || pixelCoord[1] < 0 || pixelCoord[1] > volume.getDimY() - 1 || pixelCoord[2] < 0 || pixelCoord[2] > volume.getDimZ() - 1) {
            return 0;
        }
        
        return gradient.getGradientMag();
    }
    
    /* Draw bounding box */
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

    @Override
    public void visualize(GL2 gl) {

        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);
        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);
        long startTime = System.currentTimeMillis();

        // Apply the currently selected function
        switch (function) {
	        case SLICER:
	            slicer(viewMatrix);
	            break;
	        case MIP:
	            mip(viewMatrix);
	            break;
                case COMPOSITING:
	            compositing(viewMatrix);
	            break;
                case TRANSFER_FUNCTION_2D:
                    transferFunction2D(viewMatrix);
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

        // Draw rendered image as a billboard texture
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
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
    
    /* Set renderer function */
    public void setFunction(int function) {
        this.function = function;
    }
    
    /* Enable or disable sgading option */
    public void setShading(boolean shading) {
        this.shading = shading;
    }
    
    /* Enable or disable shading option */
    public boolean getShading() {
        return this.shading;
    }
}
