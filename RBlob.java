import ij.*;

import java.util.Arrays;
import javax.vecmath.Point3d;

public class RBlob extends Point3d{
public double x,y,z;
public double mean, volume;
public double assignedCount, assignedMean, assignedVolume;
	
	public RBlob(double x,double y,double z,double mean,double volume){
		this.x = x;
		this.y = y;
		this.z = z;
		this.mean = mean;
		this.volume = volume;
		this.assignedMean = 0d;
		this.assignedVolume = 0d;
	}
	
	public RBlob(double x,double y,double z){
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public double dist(RBlob b){
	double d = Double.POSITIVE_INFINITY;
	try{
		d = ((x-b.x)*(x-b.x)) + ((y-b.y)*(y-b.y)) + ((z-b.z)*(z-b.z));
	}catch(Exception e){IJ.log(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
	return Math.sqrt(d);
	}
	
}
