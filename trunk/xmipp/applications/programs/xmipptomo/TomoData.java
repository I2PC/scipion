package xmipptomo;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;

import java.awt.Component;
import java.io.File;
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
 * extends Component because of the Java event mechanism (that is, TomoWindow can listen to this class for changes)
 */
public class TomoData extends Component {
	// name(s) of properties in the model that views may listen to 
	public static enum Properties {
		NUMBER_OF_PROJECTIONS;
	};

	
	// maybe it's better to move currentProjection to the viewer - in case different views can show different slices
	private int currentProjection=0;
	
	private static int defaultWidth=256, defaultHeight=256;
	private int width=0,height=0;
	private int numberOfProjections=0;
	
	private File file;

    private ImagePlus imp=null;
	// private java.util.List <TomoInfo> ti=Collections.emptyList();
    // private java.util.List <Float> tiltAngles=Collections.emptyList();
	// by now use a vector to store the info
	java.util.List <Float> tiltAngles=new Vector<Float> (); //imp.getStackSize());
	
	private Document tiltTextModel=null;

	// allow Views to lock waiting for the first image to load
	private Semaphore firstLoaded= new Semaphore(0);
	
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
		Float tilt=getCurrentTilt();
		if(tilt != null)
			setTiltText(tilt);
	}

	public ImagePlus getImage(){
		return imp;
	}
	
	private void setImage(ImagePlus i){
		imp=i;
	}
	
	public List<Float> getTiltAngles(){
		return tiltAngles;
	}
	
	/**
	 * @return null if tilt is undefined, tilt value otherwise
	 * 
	 */
	public Float getCurrentTilt(){
		Float t=null;
		try{
			t=getTiltAngles().get(getCurrentProjection());
		}catch (ArrayIndexOutOfBoundsException ex){
			// throw new Exception("TomoData.getCurrentTilt - undefined tilt angle");
			t=null;
		}
		return t;
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
	
	/**
	 * @return range 0..N
	 */
	public int getNumberOfProjections(){
		return numberOfProjections;
	}
	
	private void setNumberOfProjections(int n){
		if(n > 1)
			firePropertyChange(Properties.NUMBER_OF_PROJECTIONS.name(), numberOfProjections, n);
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
	
	public void import_data(String path) throws java.io.IOException, InterruptedException{	
			
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
	
	public void waitForFirstImage() throws InterruptedException{
		// Xmipp_Tomo.debug("waitForFirstImage");
		firstLoaded.acquire();
	}
	
	private void firstImageLoaded(){
		firstLoaded.release();
		// Xmipp_Tomo.debug("firstImageLoaded");
	}
	
	public void addProjection(ImageProcessor imageProcessor){
		if(getImage()==null){
			// cannot add empty stack, for example in the constructor. Better init all here			
			ImageStack stackResized=new ImageStack(imageProcessor.getWidth(), imageProcessor.getHeight());
			stackResized.addSlice(null, imageProcessor);
			setImage(new ImagePlus(getFileName(), stackResized));
		}else		
			getImage().getStack().addSlice(null, imageProcessor);
		
		setNumberOfProjections(getNumberOfProjections()+1);
		if(getNumberOfProjections() == 1)
			firstImageLoaded();
	}
}
