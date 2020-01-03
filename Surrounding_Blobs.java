import ij.*;
import ij.gui.*;
import ij.plugin.*;
import ij.process.*;
import ij.measure.*;

import java.util.Timer;
import java.util.TimerTask;
import java.util.concurrent.atomic.AtomicInteger;

import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;

import java.util.Collections;
import java.util.List;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.*;
import javax.swing.event.*;

import javax.vecmath.Point3d;

public class Surrounding_Blobs implements PlugIn, ActionListener{
private String title, Vunit;
private String path = Prefs.get("Surrounding_Blobs.path",System.getProperty("user.home"));
private JTextField C1VField,C2VField,radiusField;
private double C1minV = Prefs.get("Surrounding_Blobs.C1minV",0.01);
private double C2minV = Prefs.get("Surrounding_Blobs.C2minV",0.01);
private double radius = Prefs.get("Surrounding_Blobs.radius",0.5);	//radius to assign C1 RBlobs to C2 RBlobs
private static final String fileRegex = "^.+_[0-9]+h_.+\\.dv";
private List<Integer> timePoints;
private List<ResultsTable> C2results;
private AtomicInteger count;
private Timer timer;
private static final Font labelFont = new Font(Font.SANS_SERIF,Font.PLAIN,8);

	public void run(String _){
	try{
		JFrame gui = new JFrame("Surrounding_Blobs");
		gui.setIconImage(Toolkit.getDefaultToolkit().getImage(getClass().getResource("logo_icon.gif")));
		gui.setLayout(new BorderLayout());
		
		JPanel fieldPan = new JPanel();
		fieldPan.setLayout(new GridLayout(0,2,2,2));
		fieldPan.add(new JLabel("C1 min volume (µm³)",JLabel.RIGHT));
		C1VField = new JTextField(""+C1minV,3);
		fieldPan.add(C1VField);
		fieldPan.add(new JLabel("C2 min volume (µm³)",JLabel.RIGHT));
		C2VField = new JTextField(""+C2minV,3);
		fieldPan.add(C2VField);
		fieldPan.add(new JLabel("Assignment radius (µm)",JLabel.RIGHT));
		radiusField = new JTextField(""+radius,3);
		fieldPan.add(radiusField);
		gui.add(BorderLayout.CENTER,fieldPan);
		
		JPanel buttonPan = new JPanel();
		JButton current = new JButton("Current");
		current.addActionListener(this);
		buttonPan.add(current);
		JButton batch = new JButton("Batch");
		batch.addActionListener(this);
		buttonPan.add(batch);
		gui.add(BorderLayout.SOUTH,buttonPan);
		gui.pack();
		gui.setLocationRelativeTo(null);
		gui.setVisible(true);
	}catch(Exception e){IJ.log(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
	}
	
	//check dimensions, remove artifactual values from OMX black box processing and convert to 16-bit
	private boolean setup(ImagePlus imp){
		String title = imp.getTitle();
		//imp.setOverlay(new Overlay());
		int W = imp.getWidth();
		int H = imp.getHeight();
		int C = imp.getNChannels();
		int Z = imp.getNSlices();
		int T = imp.getNFrames();
		if(C!=2||Z==1||T>1){IJ.error(title+" dimensions are not supported");return false;}
		
		double min = 2500d;
		double max = new StackStatistics(imp).max;
		IJ.setMinAndMax(imp, min, max);
		IJ.run(imp, "16-bit", "");
		return true;
	}
	
	private void analyse(ImagePlus imp,boolean current){
	try{
		Duplicator dup = new Duplicator();
		Overlay ol = new Overlay();
		String title = imp.getTitle();
		IJ.run("Set Measurements...", "area mean standard centroid stack limit display redirect=None decimal=3");
		IJ.setBackgroundColor(0, 0, 0);
		
		ImagePlus C2 = dup.run(imp, 2, 2, 1, imp.getNSlices(), 1, 1);
		Calibration cal = imp.getCalibration();
		Vunit = " ("+cal.getUnit()+"³)";
		IJ.setAutoThreshold(C2, "Huang dark stack");
		int thresh = (int)C2.getProcessor().getMinThreshold();
		double radC2 = Math.cbrt(C2minV/((4d/3d)*Math.PI));
		double radC1 = Math.cbrt(C1minV/((4d/3d)*Math.PI));

		RBlob[] C2Blobs = new Maxima3D(C2).getBlobs(radC2,2000,thresh,C2minV);
		for(int b=0;b<C2Blobs.length;b++){
			int r = (int)Math.round(0.24d/cal.pixelWidth);
			ShapeRoi vol = new ShapeRoi(new OvalRoi((C2Blobs[b].x/cal.pixelWidth)-r,(C2Blobs[b].y/cal.pixelHeight)-r,(2*r)+1,(2*r)+1));
			int V = 0;
			double sum = 0d;
			for(int z=1;z<=imp.getNSlices();z++){
				C2.setPosition(1,z,1);
				IJ.run(C2, "Create Selection", "");
				if(C2.getRoi()==null){continue;}
				ShapeRoi tr = new ShapeRoi(C2.getRoi()).and(vol);
				if(current){
					tr.setStrokeColor(Color.MAGENTA);
					tr.setPosition(0,z,1);
					ol.add(tr);
				}
				C2.setRoi(tr);
				if(C2.getRoi()==null){continue;}
				int count = C2.getStatistics().pixelCount;
				V += count;
				sum += C2.getStatistics().mean*count;
			}
			C2Blobs[b].volume = V*cal.pixelWidth*cal.pixelHeight*cal.pixelDepth;
			C2Blobs[b].mean = sum/V;
		}
		
		ImagePlus C1 = dup.run(imp, 1, 1, 1, imp.getNSlices(), 1, 1);
		IJ.setAutoThreshold(C1, "Huang dark stack");
		thresh = (int)C1.getProcessor().getMinThreshold();
		RBlob[] C1Blobs = new Maxima3D(C1).getBlobs(radC1,2000,thresh,C1minV);
		for(int b=0;b<C1Blobs.length;b++){
			int r = (int)Math.round(0.12d/cal.pixelWidth);
			ShapeRoi vol = new ShapeRoi(new OvalRoi((C1Blobs[b].x/cal.pixelWidth)-r,(C1Blobs[b].y/cal.pixelHeight)-r,(2*r)+1,(2*r)+1));
			int V = 0;
			double sum = 0d;
			for(int z=1;z<=imp.getNSlices();z++){
				C1.setPosition(0,z,1);
				IJ.run(C1, "Create Selection", "");
				if(C1.getRoi()==null){continue;}
				ShapeRoi tr = new ShapeRoi(C1.getRoi()).and(vol);
				if(current){
					tr.setStrokeColor(Color.CYAN);
					tr.setPosition(0,z,1);
					ol.add(tr);
				}
				C1.setRoi(tr);
				if(C1.getRoi()==null){continue;}
				int count = C1.getStatistics().pixelCount;
				V += count;
				sum += C1.getStatistics().mean*count;
			}
			C1Blobs[b].volume = V*cal.pixelWidth*cal.pixelHeight*cal.pixelDepth;
			C1Blobs[b].mean = sum/V;
		}
		
		
		if(current){
			imp.setOverlay(ol);
			for(int b=0;b<C2Blobs.length;b++){
				TextRoi label = new TextRoi(C2Blobs[b].x/cal.pixelWidth,C2Blobs[b].y/cal.pixelHeight,""+b,labelFont);
				label.setStrokeColor(Color.MAGENTA);
				label.setPosition(0,(int)Math.round(C2Blobs[b].z/cal.pixelDepth),1);
				ol.add(label);
			}
		}
		
		for(int p1=0;p1<C1Blobs.length;p1++){
			double minD = Double.POSITIVE_INFINITY;
			int mindex = -1;
			for(int p2=0;p2<C2Blobs.length;p2++){
				double dist = C2Blobs[p2].dist(C1Blobs[p1]);
				if(dist<radius&&dist<minD){
					minD = dist;
					mindex = p2;
				}
			}
			if(mindex!=-1){
				C2Blobs[mindex].assignedCount++;
				C2Blobs[mindex].assignedMean += C1Blobs[p1].mean;
				C2Blobs[mindex].assignedVolume += C1Blobs[p1].volume;
			}
		}
		
		ResultsTable C2table = new ResultsTable();
		C2table.showRowNumbers(false);
		for(int p2=0;p2<C2Blobs.length;p2++){
			C2Blobs[p2].assignedMean /= C2Blobs[p2].assignedCount;
			C2Blobs[p2].assignedVolume /= C2Blobs[p2].assignedCount;
			C2table.setValue("Index",p2,p2);
			C2table.setValue("X",p2,C2Blobs[p2].x);
			C2table.setValue("Y",p2,C2Blobs[p2].y);
			C2table.setValue("Z",p2,C2Blobs[p2].z);
			C2table.setValue("Mean",p2,C2Blobs[p2].mean);
			C2table.setValue("Volume"+Vunit,p2,C2Blobs[p2].volume);
			C2table.setValue("Assigned C1 Count",p2,C2Blobs[p2].assignedCount);
			C2table.setValue("Assigned C1 Mean Intensity",p2,C2Blobs[p2].assignedMean);
			C2table.setValue("Assigned C1 Mean Volume"+Vunit,p2,C2Blobs[p2].assignedVolume);
		}
		C2table.show(title+" C2 Objects");
		
		if(!current){
			C2results.add(C2table);
			String tstr = title.substring(title.indexOf("_")+1,title.indexOf("h_"));
			int time = Integer.valueOf(tstr);
			timePoints.add(time);
			IJ.log("Finished "+title);
			count.getAndDecrement();
		}
		
	}catch(Exception e){IJ.log(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
	}
	
	private void batch(){
	try{
		JFileChooser fc = new JFileChooser(path);
		fc.setDialogTitle("Directory...");
		fc.setApproveButtonText("Select");
		fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		fc.setSelectedFile(new File(path));
		if(fc.showOpenDialog(null)==JFileChooser.APPROVE_OPTION){ 
			path = fc.getSelectedFile().getAbsolutePath();
			Prefs.set("Surrounding_Blobs.path",path);
		}
		else{
			return;
		}
		C2results = Collections.synchronizedList(new ArrayList<ResultsTable>());
		timePoints = Collections.synchronizedList(new ArrayList<Integer>());
		
		IJ.log("Surrounding Blobs : "+path);
		File inputFile = new File(path);
		boolean got = false;
		count = new AtomicInteger(0);
		if(inputFile.isDirectory()){
			File[] files = inputFile.listFiles();	
			ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
			for(int f=0;f<files.length;f++){
				if(files[f].getName().matches(fileRegex)){
					final File expDir = files[f];
					Batch job = new Batch(expDir);
					exec.submit(job);
					got = true;
					count.getAndIncrement();
				}
			}
			exec.shutdown();
			
			timer = new Timer();
			TimerTask task = new TimerTask(){	//this is a stupid way to synchronise threads, but invokeAll(ArrayList<Callable>) was blocking forever
				int n = 0;
				public void run(){
					n++;
					if(count.get()==0){
						table();
						timer.cancel();
						timer.purge();
						IJ.log(path+" done in "+IJ.d2s(n/60d,1)+" mins");
					}
				}
			};
			timer.scheduleAtFixedRate(task,1000,1000);
			
		}
		if(!got){IJ.error("No experiments found",path+" does not contain any files matching "+fileRegex);}
	}catch(Exception e){IJ.log(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}	
	}
	
	private synchronized void table(){
	try{
		ResultsTable rt = new ResultsTable();
		rt.showRowNumbers(false);
		rt.setPrecision(6);
		for(int t=0;t<=100;t++){
			double blobMean = 0d;	ArrayList<Double> blobMeanValues = new ArrayList<Double>();
			double blobVolume = 0d;	ArrayList<Double> blobVolumeValues = new ArrayList<Double>();
			double assignedCount = 0d;	ArrayList<Double> assignedCountValues = new ArrayList<Double>();
			double assignedMean = 0d;	ArrayList<Double> assignedMeanValues = new ArrayList<Double>();
			double assignedVolume = 0d;	ArrayList<Double> assignedVolumeValues = new ArrayList<Double>();
			int n = 0;
			for(int p=0;p<timePoints.size();p++){
				if(timePoints.get(p)==t){
					ResultsTable results = C2results.get(p);
					if(results.getColumnHeadings().split("\t").length!=9){throw new IllegalArgumentException("wrong number of columns in table "+results.toString());}
					for(int r=0;r<results.getCounter();r++){
						blobMean += results.getValue("Mean",r);
						blobMeanValues.add(results.getValue("Mean",r));
						blobVolume += results.getValue("Volume"+Vunit,r);
						blobVolumeValues.add(results.getValue("Volume"+Vunit,r));
						double ac = results.getValue("Assigned C1 Count",r);
						assignedCount += ac;
						assignedCountValues.add(ac);
						if(!Double.isNaN(ac)&&ac>0){
							assignedMean += results.getValue("Assigned C1 Mean Intensity",r);
							assignedMeanValues.add(results.getValue("Assigned C1 Mean Intensity",r));
							assignedVolume += results.getValue("Assigned C1 Mean Volume"+Vunit,r);
							assignedVolumeValues.add(results.getValue("Assigned C1 Mean Volume"+Vunit,r));
						}
						else{
							assignedMeanValues.add(0d);
							assignedVolumeValues.add(0d);
						}
						n++;
					}
				}
			}
			if(n>0){
				blobMean /= n;
				blobVolume /= n;
				assignedCount /= n;
				assignedMean /= n;
				assignedVolume /= n;
				double blobMeanSD = 0d;
				double blobVolumeSD = 0d;
				double assignedCountSD = 0d;
				double assignedMeanSD = 0d;
				double assignedVolumeSD = 0d;
				for(int s=0;s<n;s++){
					blobMeanSD += ( (blobMean-blobMeanValues.get(s)) * (blobMean-blobMeanValues.get(s)) );
					blobVolumeSD += ( (blobVolume-blobVolumeValues.get(s)) * (blobVolume-blobVolumeValues.get(s)) );
					assignedCountSD += ( (assignedCount-assignedCountValues.get(s)) * (assignedCount-assignedCountValues.get(s)) );
					assignedMeanSD += ( (assignedMean-assignedMeanValues.get(s)) * (assignedMean-assignedMeanValues.get(s)) );
					assignedVolumeSD += ( (assignedVolume-assignedVolumeValues.get(s)) * (assignedVolume-assignedVolumeValues.get(s)) );
				}
				blobMeanSD = Math.sqrt(blobMeanSD/n);
				blobVolumeSD = Math.sqrt(blobVolumeSD/n);
				assignedCountSD = Math.sqrt(assignedCountSD/n);
				assignedMeanSD = Math.sqrt(assignedMeanSD/n);
				assignedVolumeSD = Math.sqrt(assignedVolumeSD/n);
				
				rt.setValue("Time",rt.getCounter(),t);
				rt.setValue("C2 Blob Mean Intensity",rt.getCounter()-1,blobMean);
				rt.setValue("C2 Blob Intensity StdDev",rt.getCounter()-1,blobMeanSD);
				rt.setValue("C2 Blob Mean Volume"+Vunit,rt.getCounter()-1,blobVolume);
				rt.setValue("C2 Blob Volume StdDev",rt.getCounter()-1,blobVolumeSD);
				rt.setValue("Assigned C1 Blob Count Mean",rt.getCounter()-1,assignedCount);
				rt.setValue("Assigned C1 Blob Count StdDev",rt.getCounter()-1,assignedCountSD);
				
				rt.setValue("Assigned C1 Blob Intensity Mean",rt.getCounter()-1,assignedMean);
				rt.setValue("Assigned C1 Blob Intensity StdDev",rt.getCounter()-1,assignedMeanSD);
				rt.setValue("Assigned C1 Blob Volume Mean"+Vunit,rt.getCounter()-1,assignedVolume);
				rt.setValue("Assigned C1 Blob Volume StdDev",rt.getCounter()-1,assignedVolumeSD);
			}
		}
		rt.show("Surrounding Blobs: "+path);
	}catch(Exception e){IJ.log(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}	
	}
	
	public class Batch implements Runnable{
		private File file;
		public Batch(File file){
			try{
				this.file = file;
			}catch(Exception e){IJ.log(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}	
		}
		public void run(){
			try{
				IJ.log("Analysing "+file.getName());
				ImporterOptions BFOptions = new ImporterOptions();
				BFOptions.setColorMode(ImporterOptions.COLOR_MODE_DEFAULT);
				BFOptions.setAutoscale(true);	BFOptions.setWindowless(true);	BFOptions.setForceThumbnails(false);
				BFOptions.setConcatenate(false);	BFOptions.setGroupFiles(false);	BFOptions.setVirtual(false);
				BFOptions.setOpenAllSeries(true);
				BFOptions.setId(file.getAbsolutePath());
				ImagePlus[] images = BF.openImagePlus(BFOptions);
				for(int i=0;i<images.length;i++){
					if(setup(images[i])){
						analyse(images[i],false);
					}
					images[i].close();
				}
			}catch(Exception e){IJ.log(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}	
		}		
	}
	
	public void actionPerformed(ActionEvent ae){
		final String event = ae.getActionCommand();
		Runnable task = new Runnable(){
			public void run(){
				try{
					C1minV = Double.valueOf(C1VField.getText());
				}catch(NumberFormatException nfe){IJ.error(C1minV+" is not an number");return;}
				try{
					C2minV = Double.valueOf(C2VField.getText());
				}catch(NumberFormatException nfe){IJ.error(C2minV+" is not an number");return;}
				try{
					radius = Double.valueOf(radiusField.getText());
				}catch(NumberFormatException nfe){IJ.error(radius+" is not an number");return;}
				if(event=="Current"){
					if(WindowManager.getImageCount()==0){IJ.error("Error","What current image?");return;}
					ImagePlus imp = WindowManager.getCurrentImage();
					if(setup(imp)){
						analyse(imp,true);
					}
				}
				else if(event=="Batch"){
					batch();
				}
				Prefs.set("Surrounding_Blobs.C1minV",C1minV);
				Prefs.set("Surrounding_Blobs.C2minV",C2minV);
				Prefs.set("Surrounding_Blobs.radius",radius);
			}
		};
		SwingUtilities.invokeLater(task);
	}
	
}

