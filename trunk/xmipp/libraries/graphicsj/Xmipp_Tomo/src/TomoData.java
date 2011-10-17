/***************************************************************************
 *
 * @author: Jesus Cuenca (jcuenca@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

/**
 * - Why?
 * After choosing MVC paradigm, a Model is needed - hence TomoData
 * @extends Component because of the Java event mechanism (that is, TomoWindow can listen to this class for changes)
 */

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.StackConverter;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.File;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Semaphore;

import xmipp.MDLabel;
import xmipp.MetaData;

// TODO: remove this class when not needed
public class TomoData{
	// As a general rule, this class should not interact directly with View components: it should fire a propertyChange
	// that the View processes (modifying the required component)
	// At the same time, the Controller listens for propertyChange events from components and updates this model accordingly
	// name(s) of properties in the model that views may listen to 
	public static enum Properties {
		NUMBER_OF_PROJECTIONS,CURRENT_PROJECTION_NUMBER,CURRENT_TILT_ANGLE,CURRENT_PROJECTION_ENABLED;
	};

	private PropertyChangeSupport propertyChangeSupport;

	// PROPERTIES
	// maybe it's better to move currentProjection to the viewer - in case different views can show different slices
	// TODO: strictly speaking, the order of the projections is not specified by the sequence of the images in the
	// stack file, it is specified by the tilt angle. It's just a coincidence that the first image in the sequence is
	// the one with the smallest tilt angle.
	private int currentProjectionNumber=1;
	private int numberOfProjections=0;
	private boolean resized=false;
	
	private int width=0,height=0;
	private int originalWidth=0, originalHeight=0;
	// minimum and maximum pixel values of the projection series - for normalization
	private double min=Double.POSITIVE_INFINITY, max=Double.NEGATIVE_INFINITY;

	private File file;

    // Actual tiltseries projections are stored inside an ImagePlus
	private ImagePlus imp=null;

	// by now use a vector to store the angles - indexes can be 0..N-1
	// java.util.List <Float> tiltAngles=null;
	
	// projection enabled = preserve; projection disabled = discard
	// when using Metadata enabled field, the tilt series is not modified. Instead the selfile is.
	// java.util.List <Boolean> enabledProjections= new Vector<Boolean>();
	
	// 0.. N-1
	//List <Long> ids=new Vector<Long>();
	
	MetaData imageMetadata;

	// allow Views to wait(lock) while the first (or last) image loads
	private Semaphore firstLoaded= new Semaphore(0),lastLoaded=new Semaphore(0);
	
	public TomoData(String path){
		propertyChangeSupport = new PropertyChangeSupport(this);
		setFile(path);
	}

    public void addPropertyChangeListener(PropertyChangeListener listener) {
        propertyChangeSupport.addPropertyChangeListener(listener);
    }

    public void removePropertyChangeListener(PropertyChangeListener listener) {
        propertyChangeSupport.removePropertyChangeListener(listener);
    }

    protected void firePropertyChange(String propertyName, Object oldValue, Object newValue) {
        propertyChangeSupport.firePropertyChange(propertyName, oldValue, newValue);
    }
	
	public int getBitDepth(){
		return getImage().getBitDepth();
	}
	
	
	/**
	 * @return range 1..numberProjections
	 */
	public int getCurrentProjectionNumber() {
		return currentProjectionNumber;
	}
	
	private long getCurrentProjectionId(){
		return getStackIds() [getCurrentProjectionNumber()-1];
	}

	
	/**
	 * @param newProjectionNumber range 1..numberProjections
	 */
	public void setCurrentProjectionNumber(int newProjectionNumber) {
		
		if((newProjectionNumber<1) || (newProjectionNumber > getNumberOfProjections())){
			Logger.debug("setCurrentProjection("+newProjectionNumber+")");		
			return;
		}
		int oldProjectionNumber=currentProjectionNumber;
		this.currentProjectionNumber = newProjectionNumber;
		// possible dependencies on currentProjectionNumber should be handled by the View, not by this model
		// Unluckily, ImagePlus mixes data and presentation, so it's better to keep ImagePlus related operations here
		try{
			getImage().setSlice(getCurrentProjectionNumber());
		} catch (IllegalArgumentException ex){
			Logger.debug("setCurrentProjection - ", ex);
		}
		firePropertyChange(Properties.CURRENT_PROJECTION_NUMBER.name(), oldProjectionNumber, currentProjectionNumber);
	}

	public ImagePlus getImage(){
		return imp;
	}
	
	public long[] getStackIds(){
		long [] result=new long[getNumberOfProjections()];
		if(getMetadata() != null)
			result= getMetadata().findObjects();
		return result;
	}
	
	private void setImage(ImagePlus i){
		imp=i;
	}
	
	/* private List<Float> getTiltAngles(){
		if(tiltAngles == null){
			Vector <Float> v = new Vector<Float>();
			v.setSize(getNumberOfProjections());
			tiltAngles=v;
		}
		return tiltAngles;
	}

	private List<Boolean> getEnabledProjections(){
		return enabledProjections;
	} 
	
	public Iterator<Float> getTiltAnglesIterator(){
		return tiltAngles.iterator();
	}*/
	
	/**
	 * @return null if tilt is undefined, tilt value otherwise
	 * 
	 */
	public Double getCurrentTiltAngle(){
		Double t=null;
		try{
			t=getTiltAngle(getCurrentProjectionId());
		}catch (ArrayIndexOutOfBoundsException ex){
			t=null;
		}
		return t;
	}
	
	public Enumeration<Double> getTiltAngles(){
		int size=getStackIds().length;
		Vector<Double> tiltAngles=new Vector<Double>(size);
		for(int i=0; i< size; i++)
			tiltAngles.add(getTiltAngle(getStackIds()[i]));
		return tiltAngles.elements();
	}
	
	// Helper methods to manage tilt angles "list"
	
	/**
	 * @param i 1..N
	 */
	private Double getTiltAngle(long id){
		/*if(tiltAngles == null)
			return null;
		return getTiltAngles().get(i-1);*/
		return getMetadata().getValueDouble(MDLabel.MDL_ANGLETILT,id);
	}
	
	/**
	 * @param i 1..N
	 */
	public boolean isEnabled(long id){
		/* if(enabledProjections == null)
			return false;
		return enabledProjections.get(i-1).booleanValue();
		*/
		boolean enabled=true;
		if(getMetadata().getValueInt(MDLabel.MDL_ENABLED, id) != 1)
			enabled=false;
		return enabled;
	}
	
	public boolean isCurrentEnabled(){
		return isEnabled(getCurrentProjectionId());
	}
	
	
	/**
	 * @param i angle index - from 1 to numberOfProjections
	 * @param tilt tilt angle
	 */
	private void setTiltAngle(long id, double tilt){
		//getTiltAngles().set(i-1, new Float(tilt));
		getMetadata().setValueDouble(MDLabel.MDL_ANGLETILT, tilt, id);
	}
	
	
	/**
	 * @param i enabledProjections index - from 1 to numberOfProjections
	 * @param e true=enabled,false=disabled
	 */
	private void setEnabled(long id, boolean e){
		//getEnabledProjections().set(i-1, new Boolean(e));
		int enabled=0;
		if(e)
			enabled=1;
		getMetadata().setValueInt(MDLabel.MDL_ENABLED, enabled, id);
	}
	
	// TODO: use ids&Metadata
	public void setTiltAnglesStep(double start, double step){
		double t=start;
		for(int i=1;i<=getNumberOfProjections();i++){
			setTiltAngle(i, (float)Math.ceil(t));
			t+=step;
		}
		// setCurrentTilt(getCurrentTilt());
		firePropertyChange(Properties.CURRENT_TILT_ANGLE.name(), null, getCurrentTiltAngle());
	}
	
	public void setTiltAngles(double start, double end){
		double step= (Math.abs(end-start) / (getNumberOfProjections() - 1));
		setTiltAnglesStep(start, step);
		// no need to fire property change since setTiltAnglesStep already fires it
	}
	
	/**
	 * See if the new tilt angle series fits the number of projections 
	 * @param start
	 * @param end
	 * @param step
	 */
	public boolean checkTiltAngles(double start,double end, double step){
		int seriesCount=(int) (Math.abs(end-start) / step);
		if (seriesCount != getNumberOfProjections())
			return false;
		return true;
	}
	
	
	/**
	 * Sequential add (order matters)
	 * @param t
	 */
	public void addTiltAngle(Double t){
		// add() below needs an empty vector
		/* if(tiltAngles == null)
			tiltAngles = new Vector<Float>();
			
		getTiltAngles().add(t);
		if(tiltAngles.size() == getCurrentProjectionNumber())
			firePropertyChange(Properties.CURRENT_TILT_ANGLE.name(), null, t);*/
	}
	
	
	// TODO: there may be no tilts defined - change tilt to Float (so it can handle nulls)
	private double getInitialTilt(){
		return getCurrentTiltAngle();
	}
	
	public void addEnabled(boolean t){
		// add() below needs an empty vector
		/*if(getEnabledProjections() != null){
			getEnabledProjections().add(t);
			if(getEnabledProjections().size() == getCurrentProjectionNumber())
				firePropertyChange(Properties.CURRENT_PROJECTION_ENABLED.name(), null, t);
		}*/
	}
	
	/**
	 * @return 1..N, cyclic (next to last is 1)
	 */
	public int getNextProjection(){
		return (getCurrentProjectionNumber() % getNumberOfProjections()) + 1;
	}
	
	/**
	 * @return range 0..N
	 */
	public int getNumberOfProjections(){
		return numberOfProjections;
	}
	
	private void setNumberOfProjections(int n){
		int oldNumberOfProjections=numberOfProjections;
		numberOfProjections=n;
		if(n > 1)
			firePropertyChange(Properties.NUMBER_OF_PROJECTIONS.name(), oldNumberOfProjections, n);

	}
	
	// right now return only grayscale value as double
	// more complete method at ImagePlus.getValueAsString()
	public double getPixelValue(int x,int y){
		int[] v = getImage().getPixel(x, y);
		// 32bit images
		return Float.intBitsToFloat(v[0]);
	}
	
	
	public String getFileName(){
		return file.getName();
	}
	
	public String getFilePath(long projectionId) {
		String filename= getMetadata().getValueString(MDLabel.MDL_IMAGE,projectionId);
		// Verify filename seems correct - done in fixPath
		/* int indexAt = filename.indexOf("@");
		if(indexAt >= 0){
			// slice syntax -> should be slice_number@file_name. Remove everything before slice_number
			int i=indexAt-1;
			while(i>=0){
				if(Character.isDigit(filename.charAt(i)) == false)
					break;
				i--;
			}
			filename= filename.substring(i+1, filename.length());
		} */
		
		return filename;
	}
	
	public String getDirectory(){
		return file.getParent();
	}
	
	public String getFilePath(){
		return file.getPath();
	}

	public int getWidth() {
		return width;
	}

	public void setWidth(int width) {
		this.width = width;
	}

	public int getHeight() {
		return height;
	}

	public void setHeight(int height) {
		this.height = height;
	}
	
	public void setFile(String path){
		file=new File(path);
	}
	
	public void waitForFirstImage() throws InterruptedException{
		firstLoaded.acquire();
	}
	
	private void firstImageLoaded(){
		firstLoaded.release();
	}
	
	/**
	 * @deprecated
	 */
	public void loadCanceled(){
		firstImageLoaded();
	}
	
	public void normalize(){
		getImage().getProcessor().setMinAndMax(getMin(), getMax());
	}
	
	public void waitForLastImage() throws InterruptedException{
		lastLoaded.acquire();
	}
	
	public void lastImageLoaded(){
		lastLoaded.release();
	}
	
	/**
	 * @deprecated
	 * @param imageProcessor
	 */
	public void addProjection(ImageProcessor imageProcessor){
		if(getImage() == null){
			ImageStack stack=new ImageStack(imageProcessor.getWidth(), imageProcessor.getHeight());
			// ImagePlus does not allow to add an empty stack, hence the need to add one slice here		
			stack.addSlice(null, imageProcessor);
			setImage(new ImagePlus(getFileName(), stack));
		}else
			getImage().getStack().addSlice(null, imageProcessor);
		
		setNumberOfProjections(getNumberOfProjections()+1);
		if(getNumberOfProjections() == 1)
			firstImageLoaded();
	}
	
	public void addProjection(ImagePlus imagePlus){
		ImageProcessor imageProcessor=imagePlus.getProcessor();
		if(getImage() == null){
			ImageStack stack=new ImageStack(imageProcessor.getWidth(), imageProcessor.getHeight());
			// ImagePlus does not allow to add an empty stack, hence the need to add one slice here		
			stack.addSlice(null, imageProcessor);
			setImage(new ImagePlus(getFileName(), stack));
		}else
			getImage().getStack().addSlice(null, imageProcessor);
		
		setNumberOfProjections(getNumberOfProjections()+1);
		updateMinMax(imagePlus);
		
		if(getNumberOfProjections() == 1)
			firstImageLoaded();
	}
	
	public boolean isResized() {
		return resized;
	}
	
	/**
	 * @deprecated
	 * @return
	 */
	public boolean isLoadCanceled(){
//		return(window.getLastCommandState() == Command.State.CANCELED);
		return false;
	}


	public void setResized(boolean resized) {
		this.resized = resized;
	}
	

	public double getMin() {
		return min;
	}

	public void setMin(double min) {
		this.min = min;
	}

	public double getMax() {
		return max;
	}

	public void setMax(double max) {
		this.max = max;
	}
	
	public void updateMinMax(ImagePlus image){
		double min=image.getProcessor().getMin(), max=image.getProcessor().getMax();
		// Xmipp_Tomo.debug("updateMinMax - "+ min + "," + max);
		if(min < getMin())
			setMin(min);
		if(max > getMax())
			setMax(max);
	}
	
	public void discardCurrentProjection(){
		if(getNumberOfProjections()>1){
			setEnabled(getCurrentProjectionId(), false);
		}
		firePropertyChange(Properties.CURRENT_PROJECTION_ENABLED.name(), true, false);
	}
	
	public void enableCurrentProjection(){
		if(getNumberOfProjections()>1){
			setEnabled(getCurrentProjectionId(), true);
		}
		firePropertyChange(Properties.CURRENT_PROJECTION_ENABLED.name(), false, true);
	}
	
	/**
	 * Warning - protect your windows before calling this method
	 * @param bitDepth
	 */
	public void convertTo(int bitDepth){
		// TODO converTo - check if equation differs from ImageJ convert algorithm
		if(bitDepth == getBitDepth())
			return;
		
		switch (bitDepth) {
			case 8:
				new StackConverter(getImage()).convertToGray8();
				break;
			case 32:
				new StackConverter(getImage()).convertToGray32();
				break;
			default:
				break;
			
		}
		
			
	}
	
	public String getCurrentProjectionInfo(){
		return "Projection #" + getCurrentProjectionNumber() + " ("+ getCurrentMetadata();
	}
	public int getOriginalWidth() {
		return originalWidth;
	}

	public void setOriginalWidth(int originalWidth) {
		this.originalWidth = originalWidth;
	}

	public int getOriginalHeight() {
		return originalHeight;
	}

	public void setOriginalHeight(int originalHeight) {
		this.originalHeight = originalHeight;
	}

	private void setMetadata(MetaData md){
		imageMetadata=md;
	}
	
	public MetaData getMetadata(){
		if(imageMetadata==null)
			imageMetadata = new MetaData();
		return imageMetadata;
	}
	
	// TODO: verify possible bug: first tilt is not displayed automatically (requires scrolling). Only happens with selfiles
	/**
	 * Warning: metadata.read raises an exception when reading .sel files for which you don't have write permission
	 */
	public void readMetadata(String path){
		try {
		if(TiltSeriesIO.isSelFile(path))
			getMetadata().read(path);
		if(TiltSeriesIO.isImage(path)){
			getMetadata().read(path);
			String tltFilePath = TiltSeriesIO.getTiltFilePath(path);
			readMetadata(tltFilePath);
		}
		
		if(TiltSeriesIO.isTltFile(path))
			TiltSeriesIO.readTiltAngles(path, getMetadata());
	
		}catch (Exception ex){
			Logger.debug("readMetadata - problem with "+path+". Please check that the file exists and its permissions", ex);
		}
		firePropertyChange(Properties.CURRENT_TILT_ANGLE.name(), null, 0.0);
		firePropertyChange(Properties.CURRENT_PROJECTION_ENABLED.name(), null, false);
	}
	
	/* private Object getImageMetadata(long id,int label){
		return getMetadata().get(Long.valueOf(id)).get(Integer.valueOf(label));
	}
	
	private void setImageMetadata(long id,int label,Object value){
		Hashtable<Integer,Object> row=imageMetadata.get(id);
		if(row == null)
			row = new Hashtable<Integer,Object>();
		
		row.put(Integer.valueOf(label),value);
		imageMetadata.put(Long.valueOf(id),row);
	}
	
	public void copyMetadata(MetaData md,long id){
		ids.add(id);
		// filename
		setImageMetadata(id, MDLabel.MDL_IMAGE, md.getValueString(MDLabel.MDL_IMAGE, id));
		setImageMetadata(id, MDLabel.MDL_ANGLETILT, md.getValueDouble(MDLabel.MDL_ANGLETILT, id));

		boolean enabled=true;
		if(md.getValueInt(MDLabel.MDL_ENABLED, id) != 1)
			enabled=false;
		setImageMetadata(id, MDLabel.MDL_ENABLED, enabled);
	} */
	
	private String getCurrentMetadata(){
		long id=getCurrentProjectionId();
		return "(" + getMetadata().getValueString(MDLabel.MDL_IMAGE, id) + "). Enabled: " + getMetadata().getValueInt(MDLabel.MDL_ENABLED, id);
	}
}
