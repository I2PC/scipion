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
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Semaphore;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;

// TODO: refactor propertyChangeSupport into AbstractModel?
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
	java.util.List <Float> tiltAngles=null;
	
	// projection enabled = preserve; projection disabled = discard
	// when using Metadata enabled field, the tilt series is not modified. Instead the selfile is.
	java.util.List <Boolean> enabledProjections= new Vector<Boolean>();
	
	// This class also saves the models of the texfields of its views, one Document per textfield
	

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

	
	/**
	 * @param newProjectionNumber range 1..numberProjections
	 */
	public void setCurrentProjectionNumber(int newProjectionNumber) {
		
		if((newProjectionNumber<1) || (newProjectionNumber > getNumberOfProjections())){
			Xmipp_Tomo.debug("setCurrentProjection("+newProjectionNumber+")");		
			return;
		}
		int oldProjectionNumber=currentProjectionNumber;
		this.currentProjectionNumber = newProjectionNumber;
		// possible dependencies on currentProjectionNumber should be handled by the View, not by this model
		// Unluckily, ImagePlus mixes data and presentation, so it's better to keep ImagePlus related operations here
		try{
			getImage().setSlice(getCurrentProjectionNumber());
		} catch (IllegalArgumentException ex){
			Xmipp_Tomo.debug("setCurrentProjection - ", ex);
		}
		firePropertyChange(Properties.CURRENT_PROJECTION_NUMBER.name(), oldProjectionNumber, currentProjectionNumber);
	}

	public ImagePlus getImage(){
		return imp;
	}
	
	private void setImage(ImagePlus i){
		imp=i;
	}
	
	private List<Float> getTiltAngles(){
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
	}
	
	/**
	 * @return null if tilt is undefined, tilt value otherwise
	 * 
	 */
	public Float getCurrentTiltAngle(){
		Float t=null;
		try{
			t=getTiltAngle(getCurrentProjectionNumber());
		}catch (ArrayIndexOutOfBoundsException ex){
			t=null;
		}
		return t;
	}
	
	// Helper methods to manage tilt angles "list"
	
	/**
	 * @param i 1..N
	 */
	private Float getTiltAngle(int i){
		if(tiltAngles == null)
			return null;
		return getTiltAngles().get(i-1);
	}
	
	/**
	 * @param i 1..N
	 */
	public boolean isEnabled(int i){
		if(enabledProjections == null)
			return false;
		return enabledProjections.get(i-1).booleanValue();
	}
	
	public boolean isCurrentEnabled(){
		return isEnabled(getCurrentProjectionNumber());
	}
	
	
	/**
	 * @param i angle index - from 1 to numberOfProjections
	 * @param tilt tilt angle
	 */
	private void setTiltAngle(int i, float tilt){
		getTiltAngles().set(i-1, new Float(tilt));	
	}
	
	
	/**
	 * @param i enabledProjections index - from 1 to numberOfProjections
	 * @param e true=enabled,false=disabled
	 */
	private void setEnabled(int i, boolean e){
		getEnabledProjections().set(i-1, new Boolean(e));
	}
	
	public void setTiltAnglesStep(double start, double step){
		float t=(float)start;
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
	public void addTiltAngle(Float t){
		// add() below needs an empty vector
		if(tiltAngles == null)
			tiltAngles = new Vector<Float>();
			
		getTiltAngles().add(t);
		if(tiltAngles.size() == getCurrentProjectionNumber())
			firePropertyChange(Properties.CURRENT_TILT_ANGLE.name(), null, t);
	}
	
	public void deleteTiltAngle(int i){
		if(getTiltAngles() != null)
			if( (i < getTiltAngles().size()) && (i> 1))
				getTiltAngles().remove(i-1);
	}
	
	
	// TODO: there may be no tilts defined - change tilt to Float (so it can handle nulls)
	private float getInitialTilt(){
		return getCurrentTiltAngle();
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
	
	// TODO
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
	
	public void addProjection(ImageProcessor imageProcessor){
		if(getImage() == null)
			createImageStack(imageProcessor);
		getImage().getStack().addSlice(null, imageProcessor);
		
		// TODO: addProjection - should be the value read from metadata
		getEnabledProjections().add(new Boolean(true));
		
		setNumberOfProjections(getNumberOfProjections()+1);
		if(getNumberOfProjections() == 1)
			firstImageLoaded();
	}
	
	private void createImageStack(ImageProcessor imageProcessor){
		// cannot add empty stack, for example in the constructor. better init all here			
		ImageStack stack=new ImageStack(imageProcessor.getWidth(), imageProcessor.getHeight());
		stack.addSlice(null, imageProcessor);
		setImage(new ImagePlus(getFileName(), stack));
		// getImage().setWindow(window);
	}


	public boolean isResized() {
		return resized;
	}
	
	// TODO: isLoadCanceled - manage with thread.cancel
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
			setEnabled(getCurrentProjectionNumber(), false);
		}
		firePropertyChange(Properties.CURRENT_PROJECTION_ENABLED.name(), true, false);
	}
	
	public void enableCurrentProjection(){
		if(getNumberOfProjections()>1){
			setEnabled(getCurrentProjectionNumber(), true);
		}
		firePropertyChange(Properties.CURRENT_PROJECTION_ENABLED.name(), false, true);
	}
	
	/**
	 * Warning - protect your windows before calling this method
	 * @param bitDepth
	 */
	public void convertTo(int bitDepth){
		// TODO converTo - check if Coss equation fits better than ImageJ convert algorithm
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
		return "Projection #" + getCurrentProjectionNumber() + ", enabled:" + isCurrentEnabled();
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

}
