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
    
    // Shading option
    private boolean shading;
    
    // Phong parameters 
    private double kAmbient = 0.1;
    private double kDiffuse = 0.7;
    private double kSpecular = 0.2;
    private double kAlpha = 10;
    
    private double[] phongParams = new double[] {
            this.kAmbient,
            this.kDiffuse,
            this.kSpecular,
            this.kAlpha
    };
    
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
    
    // This function linearly interpolates the value p0 and p1 given the factor between 0 and 1 
    private float linearInterpolate(float p0, float p1, float factor) {
        
    	return p0 * (1 - factor) + p1 * factor; 
    }
    
    // This function computes the intensity value for the point with coordinates coord, using trilinear interpolation
    private short getVoxelLinearInterpolated(double[] coord) {
        
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
        the factors (between 0 and 1) needed for interpolate function are (coord-x0)/(x1-x0) = coord-x0 because x1-x0=1.*/
	  
	/* Performing trilinear interpolation */
        
        // First we interpolate along x (note: coding used for values is cyz, so c10 means y=y0+1, z=z0)
        // Bottom cube face
	float c00 = linearInterpolate(volume.getVoxel(x0, y0, z0), volume.getVoxel(x0 + 1, y0, z0), xd);
	float c10 = linearInterpolate(volume.getVoxel(x0, y0 + 1, z0), volume.getVoxel(x0 + 1, y0 + 1, z0), xd);
	// Top cube face
	float c01 = linearInterpolate(volume.getVoxel(x0, y0, z0 + 1), volume.getVoxel(x0 + 1, y0, z0 + 1), xd);
	float c11 = linearInterpolate(volume.getVoxel(x0, y0 + 1,z0 + 1), volume.getVoxel(x0 + 1, y0 + 1, z0 + 1), xd); 
	    
	// Interpolate along y (note: v0/v1 refers to whether z=z0 or z=z0+1)
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
                int val = this.interactiveMode ? getVoxel(pixelCoord) : getVoxelLinearInterpolated(pixelCoord);
                
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
        
        // Array to store the entry and exit point of the ray
        int[] range = new int[2];
        int entryPoint;
        int exitPoint;
        
        // Ray computation for each pixel
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {

                // Initialize max intensity
                int maxIntensity = 0;
                
                // Compute the optimal entry and exit point of the ray for the current pixel
                range = computeRayRange(i, j, viewVec, uVec, vVec, imageCenter);
                entryPoint = range[0];
                exitPoint = range[1];
                
                // Steps along the ray to compute the mazimum intensity projection for the current pixel
                for (int k = entryPoint; k < exitPoint; k+=step) {
             
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * k + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * k + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * k + volumeCenter[2];

                    // Speed up rendering when moving volume
                    int val = this.interactiveMode ? getVoxel(pixelCoord) : getVoxelLinearInterpolated(pixelCoord);
                    
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
    
    /* Compute the optimal entry and exit point of the ray */
    private int[] computeRayRange(int i, int j, double[] viewVec, double[] uVec, double[] vVec, int imageCenter) {
        
        // Initialize entry and exit point of the ray
        int start = Integer.MIN_VALUE;
        int end = Integer.MAX_VALUE;
        
        // Volume dimensions
        int[][] volumeBoundaries = new int[3][2];
        volumeBoundaries[0][0] = - volume.getDimX() / 2;
        volumeBoundaries[0][1] = volume.getDimX() / 2;
        volumeBoundaries[1][0] = - volume.getDimY() / 2;
        volumeBoundaries[1][1] = volume.getDimY() / 2;
        volumeBoundaries[2][0] = - volume.getDimZ() / 2;
        volumeBoundaries[2][1] = volume.getDimZ() / 2;

        // Origin 
        double[] origin = new double[3];
        
        // Range [entryPoint, exitPoint]
        int[] range = new int[2];

        for (int c = 0; c < 3; c++) {

            // Set origin
            origin[c] = (i - imageCenter) * uVec[c] + (j - imageCenter) * vVec[c];

            // If the origin is not between the boundaries then return range = [0, 0]
            if (viewVec[c] == 0) {
                
                if ((origin[c] < volumeBoundaries[c][0] || origin[c] > volumeBoundaries[c][1])) {
                    
                    return new int[]{0, 0};
                }
            } else {

                // Otherwise for each dimension find low/high: range[x,y,z][low, high] : range[0, 1, 2][0, 1]
                range[0] = (int) ((volumeBoundaries[c][0] - origin[c]) / viewVec[c]);
                range[1] = (int) ((volumeBoundaries[c][1] - origin[c]) / viewVec[c]);

                // If low > high, swap
                if (range[0] > range[1]) {
                    
                    int tmp = range[0];
                    range[0] = range[1];
                    range[1] = tmp;
                }
                // If low > start , start = low
                if (range[0] > start) {
                    
                    start = range[0];
                }
                // If high < end , end = high
                if (range[1] < end) {
                    
                    end = range[1];
                }
                
                if (start > end || end < 0) {
                    
                    return new int[]{0, 0};
                }
            }
        }
        
        // Return entry and exit point of the ray
        return new int[]{start, end};
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

        short[][] sumIntensity = new short[image.getHeight()][image.getWidth()];
        
        // Array to store the entry and exit point of the ray
        int[] range = new int[2];
        int entryPoint;
        int exitPoint;

         // Step size of samples akong the viewing ray
        step = this.interactiveMode ? INTERACTIVE_MODE_STEP : NON_INTERACTIVE_MODE_STEP;

        // Ray computation for each pixel
        for (int j = 0; j < image.getHeight(); j ++) {
            for (int i = 0; i < image.getWidth(); i ++) {
                
                // Initialize sumIntensity
                sumIntensity[i][j] = 0;
                
                // Compute the optimal entry and exit point of the ray for the current pixel
                range = computeRayRange(i, j, viewVec, uVec, vVec, imageCenter);
                entryPoint = range[0];
                exitPoint = range[1];
        
                TFColor color = new TFColor(0, 0, 0, 0);
                
                /* Steps along the ray for the current pixel */
                for (int k = entryPoint; k < exitPoint; k+=step) {
                    
                    // Get calculate new volumeCenter
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * k + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * k + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * k + volumeCenter[2];

                    // Speed up rendering when moving volume
                    int val = this.interactiveMode ? getVoxel(pixelCoord) : getVoxelLinearInterpolated(pixelCoord);

                    voxelColor = tFunc.getColor(val);

                    color.r = voxelColor.r * voxelColor.a + (1 - voxelColor.a) * color.r;
                    color.g = voxelColor.g * voxelColor.a + (1 - voxelColor.a) * color.g;
                    color.b = voxelColor.b * voxelColor.a + (1 - voxelColor.a) * color.b;

                    color.a = (1 - voxelColor.a) * color.a;

                }
                color.a = 1 - color.a;
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = color.a <= 1.0 ? (int) Math.floor(color.a * 255) : 255;
                int c_red = color.r <= 1.0 ? (int) Math.floor(color.r * 255) : 255;
                int c_green = color.g <= 1.0 ? (int) Math.floor(color.g * 255) : 255;
                int c_blue = color.b <= 1.0 ? (int) Math.floor(color.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                image.setRGB(i, j, pixelColor);
            }
        }
    }
    
    
    /* 2D transfer function */
    void transferFunction2D(double[] viewMatrix) {

        // Clear image
        this.clearImage();
        
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        double[] reverseView = new double[3];
        VectorMath.setVector(reverseView, -viewMatrix[2], -viewMatrix[6], -viewMatrix[10]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        TFColor voxelColor = new TFColor();

        short[][] sumIntensity = new short[image.getHeight()][image.getWidth()];
        
        int[] range = new int[2];
        int entryPoint;
        int exitPoint;
        double[] rgb;

        // Step size of samples akong the viewing ray
        step = this.interactiveMode ? INTERACTIVE_MODE_STEP : NON_INTERACTIVE_MODE_STEP;
        
        // 2D transfer function computation for each pixel
        for (int j = 0; j < image.getHeight(); j ++) {
            for (int i = 0; i < image.getWidth(); i ++) {
                
                // Initialize sumIntensity
                sumIntensity[i][j] = 0;
                
                // Compute the optimal entry and exit point of the ray for the current pixel
                range = computeRayRange(i,j, viewVec, uVec, vVec, imageCenter);
                entryPoint = range[0];
                exitPoint = range[1];

                TFColor compositingColor = new TFColor(0, 0, 0, 0);
                VoxelGradient gradient;
                
                // Steps alon the ray for the current pixel
                for (int k = exitPoint; k > entryPoint; k-=step) {

                    // Get calculate new volumeCenter
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * (k) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * (k) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * (k) + volumeCenter[2];

                    // Speed up rendering when moving volume
                    int val = this.interactiveMode ? getVoxel(pixelCoord) : getVoxelLinearInterpolated(pixelCoord);

                    voxelColor = tfEditor2D.triangleWidget.color.clone();

                    gradient = gradients.getTriLinearGradient((float) pixelCoord[0], (float) pixelCoord[1], (float) pixelCoord[2]);
                    double dotProduct = VectorMath.dotproduct(viewVec, gradient.normalisedVector());
                    double opacity = computeOpacity(val, gradient);
                    
                    if (this.shading) { // Shading option
                        
                        rgb = new double[]{0, 0, 0};
                        
                        if (gradient.mag > 0 && dotProduct > 0 && opacity > 0) {
                            
                            double[] compRGB = new double[]{voxelColor.r, voxelColor.g, voxelColor.b};
                            
                            for (int z = 0; z < 3; z++) {
                                
                                rgb[z] = phongParams[0] + compRGB[z] * phongParams[1] * dotProduct + phongParams[2] * Math.pow(dotProduct, phongParams[3]);
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
            }
        }
    }

    /* Compute the opacity for the 2D transfer function */
    private double computeOpacity(double intensity, VoxelGradient voxelGradient) {
        
        TFColor color = tfEditor2D.triangleWidget.color;
        short baseIntensity = tfEditor2D.triangleWidget.baseIntensity;
        double radius = tfEditor2D.triangleWidget.radius;

        float gradientMagnitude = voxelGradient.mag;

        if (voxelGradient.mag < tfEditor2D.triangleWidget.minGradient || voxelGradient.mag > tfEditor2D.triangleWidget.maxGradient) {
            
            return 0;
        }

        if (gradientMagnitude == 0 && intensity == baseIntensity) {
            
            return color.a;
        } else if (gradientMagnitude > 0 && (intensity - radius * gradientMagnitude <= baseIntensity) && (baseIntensity <= intensity + radius * gradientMagnitude)) {
            
            return color.a * (1 - (1 / radius) * Math.abs((baseIntensity - intensity) / (gradientMagnitude)));
        }
        
        return 0.0;
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
}
