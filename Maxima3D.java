import ij.*;
import ij.gui.*;
import ij.plugin.*;
import ij.process.*;
import ij.measure.*;

import java.util.ArrayList;
import java.util.Arrays;

import javax.vecmath.Point3d;

public class Maxima3D{
private int W,H,C,Z,T;
private double pixelW,pixelD,zRatio;
private int[][][] vox;

public Maxima3D(ImagePlus imp){
try{
	this.W = imp.getWidth();
	this.H = imp.getHeight();
	this.C = imp.getNChannels();
	this.Z = imp.getNSlices();
	this.T = imp.getNFrames();
	if(C>1||Z==1||T>1){throw new IllegalArgumentException("3D z-stack required");}
	this.vox = new int[W][H][Z];
		IJ.run(imp, "Select None", "");
		for(int z=0;z<Z;z++){
			imp.setSlice(z+1);
			vox[z]= imp.getProcessor().getIntArray();
		}
	Calibration cal = imp.getCalibration();
	this.pixelW = cal.pixelWidth;
	this.pixelD = cal.pixelDepth;
	this.zRatio = pixelD/pixelW + 1d;
}catch(Exception e){IJ.log(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
}

public RBlob[] getBlobs(double radius, int tolerance, int threshold, double minV){
	RBlob[] blobs = new RBlob[0];
try{
	Point3d[] points = this.findMaxima(radius, tolerance, threshold, minV);
	blobs = new RBlob[points.length];
	for(int b=0;b<points.length;b++){
		blobs[b] = new RBlob(points[b].x,points[b].y,points[b].z);
	}
}catch(Exception e){IJ.log(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
	return blobs;
}

public Point3d[] findMaxima(double radius, int tolerance, int threshold, double minV){
	Point3d[] points = new Point3d[0];
	try{
		ArrayList<Point3d> pList = new ArrayList<Point3d>();
		int minCount = (int)Math.round(minV/pixelW/pixelW/pixelD);
		//radius = calibrated sampling radius
		double radius2 = radius * radius;	//square radius for comparison to avoid sqrt in ED
		int rXY = (int)Math.ceil(radius/pixelW);	//3D sampling box radii
		int rZ = (int)Math.ceil(radius/pixelD);
		rZ = (int)Math.round(rZ*zRatio);
		for(int z=rZ;z<Z-rZ;z++){
			IJ.showStatus("Calculating 3D maxima...");
			IJ.showProgress(z,Z-rZ);
			for(int x=rXY;x<W-rXY;x++){
				centreLoop:
				for(int y=rXY;y<H-rXY;y++){
					//**********3D patch
					int centre = vox[z][x][y];
					int count = 0;
					if(centre<threshold){
						continue centreLoop;
					}
					for(int pz=-rZ;pz<=rZ;pz++){
						for(int px=-rXY;px<=rXY;px++){
							for(int py=-rXY;py<=rXY;py++){
								if(px==0&&py==0&&pz==0){continue;}
								//prolate spheroid radius
								double radD2 = ((pz*pz*pixelD*pixelD)/zRatio) + (py*py*pixelW*pixelW) + (px*px*pixelW*pixelW) ;
								if(radD2<=radius2){
									if(vox[z+pz][x+px][y+py]>threshold){
										count++;
									}
									if(centre+tolerance<vox[z+pz][x+px][y+py]){
										continue centreLoop;
									}
								}
							}
						}
					}
					if(count>=minCount){
						pList.add(new Point3d(x*pixelW,y*pixelW,(z+1)*pixelD));	//if no continue criteria were met, this is a maximum
						vox[z][x][y] = Integer.MAX_VALUE;		//stop any other maxima being called within radius
					}
					//*********3D patch
				}
			}
		}
		IJ.showStatus("");
		
		points = new Point3d[pList.size()];
		for(int i=0;i<pList.size();i++){
			points[i] = pList.get(i);
		}
	}catch(Exception e){IJ.log(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
return points;
}

}
