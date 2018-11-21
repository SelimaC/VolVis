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

    private int function;
    
    final static int INTERACTIVE_MODE_STEP = 4;
    final static int INTERACTIVE_MODE_GRANULARITY = 4;
    final static int NON_INTERACTIVE_MODE_STEP = 1;
    final static int NON_INTERACTIVE_MODE_GRANULARITY = 1;

    // Ray parameters
    int granularity = NON_INTERACTIVE_MODE_GRANULARITY;
    int step = NON_INTERACTIVE_MODE_STEP;
    
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
        //tFunc.setTestFunc();
        
        
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
     

    private short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return volume.getVoxel(x, y, z);
    }
    
    public double getTriLinearInterpolatedVoxel(double[] pixelCoord)
    {
        double x = pixelCoord[0];
        int xhigh = (int) Math.ceil(x);
        int xlow = (int)Math.floor(x);
        double y = pixelCoord[1];
        int yhigh = (int)Math.ceil(y);
        int ylow = (int)Math.floor(y);
        double z = pixelCoord[2];
        int zhigh = (int)Math.ceil(z);
        int zlow = (int)Math.floor(z);
        if (xlow < 0 || xlow >= volume.getDimX() || ylow < 0 || ylow >= volume.getDimY()
                   || zlow < 0 || zlow >= volume.getDimZ()) {
               return 0;
        }
        if (xhigh < 0 || xhigh >= volume.getDimX() || yhigh < 0 || yhigh >= volume.getDimY()
                   || zhigh < 0 || zhigh >= volume.getDimZ()) {
               return 0;
        }

        //Pepare the volxe cube coordinate, x0 to x3 are the bottom coordinate, others are upper coordinate.
        int[] x0 = {xlow, ylow, zlow};
        int[] x1 = {xhigh, ylow, zlow};
        int[] x2 = {xhigh, yhigh, zlow};
        int[] x3 = {xhigh, yhigh, zlow};
     
        int[] x4 = {xlow,ylow,zhigh};
        int[] x5 = {xhigh,ylow,zhigh};
        int[] x6 = {xlow,yhigh,zhigh};
        int[] x7 = {xhigh, yhigh,zhigh};
        //get v value of each point
        short sx0 = volume.getVoxelbycoordinate(x0);
        short sx1 = volume.getVoxelbycoordinate(x1);
        short sx2 = volume.getVoxelbycoordinate(x2);
        short sx3 = volume.getVoxelbycoordinate(x3);
        short sx4 = volume.getVoxelbycoordinate(x4);
        short sx5 = volume.getVoxelbycoordinate(x5);
        short sx6 = volume.getVoxelbycoordinate(x6);
        short sx7 = volume.getVoxelbycoordinate(x7);
        //calculate alpha, beta, and gamma in the formula
        // length(x-x0)/length(x1-x0)
        double[] vxxo = {x-xhigh, y-ylow, z-zlow};
        double[] vx1x0 = {xhigh-xlow, ylow-ylow, zlow-zlow};
        double alpha = VectorMath.length(vxxo)/VectorMath.length(vx1x0);
        //beta = (x-x1)/(x3-x1)
        double[] vxx1 = {x-xhigh, y-ylow, z-zlow};
        double[] vx3x1 = {xhigh -xhigh, yhigh - ylow, zlow-zlow};
        double beta = VectorMath.length(vxx1)/VectorMath.length(vx3x1);
        //gamma = (x-x3/x7-x3)
        double[] vxx3 = {x-xhigh, y-yhigh, z-zlow};
        double[] vx7x3 = {xhigh - xhigh, yhigh - yhigh, zhigh - zlow};
        double gamma = VectorMath.length(vxx3)/VectorMath.length(vx7x3);
        //calculate the scalar value
        double St = (1-alpha)*(1-beta)*(1-gamma)*sx0 + alpha*(1-beta)*(1-gamma)*sx1
          +(1-alpha)*beta*(1-gamma)*sx2+ alpha*beta*(1-gamma)*sx3
          +(1-alpha)*(1-beta)*gamma*sx4+alpha*(1-beta)*gamma*sx5
          +(1-alpha)*beta*gamma*sx6+ alpha*beta*gamma*sx7 ;
        return St;
    }

    void clearImage(){
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
    }

    void slicer(double[] viewMatrix) {

        // clear image
        this.clearImage();

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

        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                int val = getVoxel(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
                
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


    void mip(double[] viewMatrix) {

        // clear image
        this.clearImage();

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

        
        int[] kRange = new int[2];
        for (int j = 0; j < image.getHeight(); j ++) {
            for (int i = 0; i < image.getWidth(); i ++) {

                int maxIntensity = 0;

                kRange = optimalDepth(imageCenter, viewVec, uVec, vVec, i, j);
                for (int k = kRange[0]; k < kRange[1]; k++) {
                    // Get calculate new volumeCenter
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * (k) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * (k) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * (k) + volumeCenter[2];

                    int val = (short) getTriLinearInterpolatedVoxel(pixelCoord);
                    if (val > maxIntensity) {
                        maxIntensity = val;
                    }
                }

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = maxIntensity / max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = maxIntensity > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque

                long pixelColor = this.pixelColor(voxelColor);
                image.setRGB(i, j, (int) pixelColor);

            }
        }

    }
    
    private long pixelColor(TFColor v) {
        int c_alpha = v.a <= 1.0 ? (int) Math.floor(v.a * 255) : 255;
        int c_red = v.r <= 1.0 ? (int) Math.floor(v.r * 255) : 255;
        int c_green = v.g <= 1.0 ? (int) Math.floor(v.g * 255) : 255;
        int c_blue = v.b <= 1.0 ? (int) Math.floor(v.b * 255) : 255;
        return this.binaryColor(new long[]{c_alpha, c_red, c_green, c_blue});
    }
    
    private long binaryColor(long[] rgba) {
        return (rgba[0] << 24) | (rgba[1] << 16) | (rgba[2] << 8) | rgba[3];
    }
    
    private int[] optimalDepth(int imageCenter, double[] viewVec, double[] uVec, double[] vVec, int i, int j) {
        int kStart = Integer.MIN_VALUE;
        int kEnd = Integer.MAX_VALUE;
        double[] pixelCoord = new double[3];

        //if (!this.planeIntersection) {
          //  return new int[]{0, volume.getDiagonalDepth()};
        //}
        //set vloume dimensions volume[x,y,z][low,high] = vloume[0,1,2][0,1]
        int[][] volumeBoundary = new int[3][2];
        volumeBoundary[0][0] = -volume.getDimX() / 2;
        volumeBoundary[0][1] = volume.getDimX() / 2;
        volumeBoundary[1][0] = -volume.getDimY() / 2;
        volumeBoundary[1][1] = volume.getDimY() / 2;
        volumeBoundary[2][0] = -volume.getDimZ() / 2;
        volumeBoundary[2][1] = volume.getDimZ() / 2;

        //origin[x,y,z]
        double[] origin = new double[3];

        int[] kList = new int[2];

        for (int l = 0; l < 3; l++) {

            //set origin  origin[x,y,z] = origin[0,1,2]
            origin[l] = (i - imageCenter) * uVec[l] + (j - imageCenter) * vVec[l];

            //if origin not between the slabs then return
            if (viewVec[l] == 0) {
                if ((origin[l] < volumeBoundary[l][0] || origin[l] > volumeBoundary[l][1])) {
                    return new int[]{0, 0};
                }
            } else {

                //for each dimension find kLow/kHigh . k[x,y,z][low,high] : k[0,1,2][0,1]
                kList[0] = (int) ((volumeBoundary[l][0] - origin[l]) / viewVec[l]);
                kList[1] = (int) ((volumeBoundary[l][1] - origin[l]) / viewVec[l]);

                // if kLow > kHigh, swap
                if (kList[0] > kList[1]) {
                    int tmp = kList[0];
                    kList[0] = kList[1];
                    kList[1] = tmp;
                }
                //if kLow > kStart , kStart = kLow
                if (kList[0] > kStart) {
                    kStart = kList[0];
                }
                //if kHigh < kEnd , kEnd = kHigh
                if (kList[1] < kEnd) {
                    kEnd = kList[1];
                }
                //
                if (kStart > kEnd || kEnd < 0) {
                    return new int[]{0, 0};
                }
            }
        }
        return new int[]{kStart, kEnd};

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
	            mip(viewMatrix);
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
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
    
    public void setFunction(int function) {
        this.function = function;
    }
}
