/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

import java.io.File;
import java.io.IOException;

/**
 *
 * @author michel
 */
public class Volume {
    
    public Volume(int xd, int yd, int zd) {
        data = new short[xd*yd*zd];
        dimX = xd;
        dimY = yd;
        dimZ = zd;
    }
    
    public Volume(File file) {
        
        try {
            VolumeIO reader = new VolumeIO(file);
            dimX = reader.getXDim();
            dimY = reader.getYDim();
            dimZ = reader.getZDim();
            data = reader.getData().clone();
            computeHistogram();
        } catch (IOException ex) {
            System.out.println("IO exception");
        }
        
    }   
    
    public short getVoxel(int x, int y, int z) {
        return data[x + dimX*(y + dimY * z)];
    }
    
    public void setVoxel(int x, int y, int z, short value) {
        data[x + dimX*(y + dimY*z)] = value;
    }

    public void setVoxel(int i, short value) {
        data[i] = value;
    }
    
    public short getVoxel(int i) {
        return data[i];
    }
    
    //This function linearly interpolates the value g0 and g1 given the factor (t) between [0,1] 
    //the interpolated value is returned.
    private float interpolate(float g0, float g1, float factor) {
    	float result= g0 * (1-factor)+g1*factor;
    	return result; 
    }
    
    // This function computes the intensity value for the point with coordinates coord, using trilinear interpolation
	public short getVoxelLinearInterpolate(double[] coord) {
	    // making sure that point coordinates are within the limits.
	    if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
	            || coord[2] < 0 || coord[2] > (dimZ-2)) {
	        return 0;
	    }
	    /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
	    
	    //Identify the cell which contains the point with coordinates coord
	    int x0 = (int) Math.floor(coord[0]); 
	    int y0 = (int) Math.floor(coord[1]);
	    int z0 = (int) Math.floor(coord[2]);
	    
	    // since we assume the distances between neighbouting voxels are 1 in all directions: 
	    // the factors (between 0 and 1) needed for interpolate function are (coord - x0)/(x1-x0) = coord - x0 because x1-x0=1.
	    float xdif = (float) (coord[0]-x0);
	    float ydif = (float) (coord[1]-y0);
	    float zdif = (float) (coord[2]-z0);
	    
	    // Performing trilinear interpolation
	    // interpolate along x (note: coding used for values is vyz, so v10 means y=y0+1, z=z0)
	    // bottom cube face
	    float v00=interpolate(getVoxel(x0,y0,z0),getVoxel(x0 + 1,y0,z0), xdif);
	    float v10=interpolate(getVoxel(x0,y0 + 1,z0),getVoxel(x0 + 1,y0 + 1,z0), xdif); 
	    
	    //top cube face
	    float v01=interpolate(getVoxel(x0,y0,z0 + 1),getVoxel(x0 + 1,y0,z0 + 1), xdif);
	    float v11=interpolate(getVoxel(x0,y0 + 1,z0 + 1),getVoxel(x0 + 1,y0 + 1,z0 + 1), xdif); 
	    
	    //interpolate along y (note: v0/v1 refers to whether z=z0 or z=z0+1)
	    float v0 =interpolate(v00,v10,ydif);
	    float v1 =interpolate(v01,v11,ydif);
	    
	    //interpolate along z
	    float vInt = interpolate(v0,v1,zdif);
	    
	    return (short)vInt; 
	}
    
    public int getDimX() {
        return dimX;
    }
    
    public int getDimY() {
        return dimY;
    }
    
    public int getDimZ() {
        return dimZ;
    }

    public short getMinimum() {
        short minimum = data[0];
        for (int i=0; i<data.length; i++) {
            minimum = data[i] < minimum ? data[i] : minimum;
        }
        return minimum;
    }

    public short getMaximum() {
        short maximum = data[0];
        for (int i=0; i<data.length; i++) {
            maximum = data[i] > maximum ? data[i] : maximum;
        }
        return maximum;
    }
 
    public int[] getHistogram() {
        return histogram;
    }
    
    private void computeHistogram() {
        histogram = new int[getMaximum() + 1];
        for (int i=0; i<data.length; i++) {
            histogram[data[i]]++;
        }
    }
    
    private int dimX, dimY, dimZ;
    private short[] data;
    private int[] histogram;
}
