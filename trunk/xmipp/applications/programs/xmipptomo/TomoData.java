package xmipptomo;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.LookUpTable;
import ij.io.FileInfo;
import ij.io.FileOpener;
import ij.io.ImageReader;
import ij.io.OpenDialog;
import ij.process.ImageProcessor;

import java.awt.image.ColorModel;
import java.awt.image.IndexColorModel;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.util.Collections;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Semaphore;

import javax.swing.text.BadLocationException;
import javax.swing.text.Document;

/**
 * XmippTomo data model: tilt series with its projections
 */

/**
 * @author jcuenca
 * 
 */

/**
 * @author jcuenca
 *
 */
public class TomoData {
	
	// maybe it's better to move currentProjection to the viewer - in case different views can show different slices
	private int currentProjection=0;
	
	private static int defaultWidth=256, defaultHeight=256;
	private int width,height;
	private int numberOfProjections;
	
	private File file;

    private ImagePlus imp=null;
	// private java.util.List <TomoInfo> ti=Collections.emptyList();
    // private java.util.List <Float> tiltAngles=Collections.emptyList();
	// by now use a vector to store the info
	java.util.List <Float> tiltAngles=new Vector<Float> (); //imp.getStackSize());
	
	private Document tiltTextModel=null;

	// allow Views to wait for the first image to load
	private Semaphore firstLoaded= new Semaphore(1);
	
	public void setTiltModel(Document model){
		if(tiltTextModel==null){
			// initialize model with default value
			try{
				model.insertString(0, String.valueOf(getInitialTilt()), null);
			}catch (Exception ex){}
		}
		
		tiltTextModel=model;
	}
	
	public Document getTiltModel(){
		return tiltTextModel;
	}
	
	
	/**
	 * @return range 0..numberProjections-1
	 */
	public int getCurrentProjection() {
		return currentProjection;
	}

	
	/**
	 * @param currentProjection range 0..numberProjections-1
	 */
	public void setCurrentProjection(int currentProjection) {
		this.currentProjection = currentProjection;
		// update all things depending on current slice
		getImage().setSlice(getCurrentProjection());
		setTiltText(getCurrentTilt());
	}

	public ImagePlus getImage(){
		return imp;
	}
	
	public void setImage(ImagePlus i){
		imp=i;
	}
	
	public List<Float> getTiltAngles(){
		return tiltAngles;
	}
	
	public float getCurrentTilt(){
		return getTiltAngles().get(getCurrentProjection()).floatValue();
	}
	
	public void setCurrentTilt(float t){
		getTiltAngles().set(getCurrentProjection(), new Float(t));
	}
	
	// method for loading the angles in sequence (it assumes an implicit order) 
	public void addTiltAngle(Float t){
		tiltAngles.add(t);
	}
	
	// there may be no tilts defined - change tilt to Float (so it can handle nulls)
	private float getInitialTilt(){
		return getCurrentTilt();
	}
	
	public int getNumberOfProjections(){
		return numberOfProjections;
	}
	
	public void setNumberOfProjections(int n){
		numberOfProjections=n;
	}
	
	// right now return only grayscale value as double
	// more complete method at ImagePlus.getValueAsString()
	public double getPixelValue(int x,int y){
		int[] v = getImage().getPixel(x, y);
		// 32bit images
		return Float.intBitsToFloat(v[0]);
		//return getImage().getCalibration().getCValue(v[0]);
	}
	
	public void import_data(String path) throws java.io.IOException{	
			
			new TiltSeriesOpener().read(path,this);

	}
	
	
	private void setDocumentText(Document d,String t){
		try{
			d.remove(0,d.getLength());
			d.insertString(0, t, null);
		}catch (BadLocationException e){}
	}
	
	private void setTiltText(float tilt){
		setDocumentText(getTiltModel(),String.valueOf(tilt));	
	}

	/**
	 * @return the defaultWidth
	 */
	public int getDefaultWidth() {
		return defaultWidth;
	}

	/**
	 * @return the defaultHeight
	 */
	public int getDefaultHeight() {
		return defaultHeight;
	}
	
	public String getFileName(){
		return file.getName();
	}

	/**
	 * @return the width
	 */
	public int getWidth() {
		return width;
	}

	/**
	 * @param width the width to set
	 */
	public void setWidth(int width) {
		this.width = width;
	}

	/**
	 * @return the height
	 */
	public int getHeight() {
		return height;
	}

	/**
	 * @param height the height to set
	 */
	public void setHeight(int height) {
		this.height = height;
	}
	
	public void setFile(String path){
		file=new File(path);
	}
	
	public synchronized void waitForFirstImage() throws InterruptedException{
		firstLoaded.acquire();
	}
	
	private synchronized void firstImageLoaded(){
		firstLoaded.release();
	}
}
