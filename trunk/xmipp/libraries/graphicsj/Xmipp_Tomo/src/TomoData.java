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

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;

import java.awt.Component;
import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Semaphore;

import javax.swing.text.BadLocationException;
import javax.swing.text.Document;


/**
 * Data model, following the MVC paradigm
 * @author jcuenca
 * extends Component because of the Java event mechanism (that is, TomoWindow can listen to this class for changes)
 */
public class TomoData extends Component {
	// name(s) of properties in the model that views may listen to 
	public static enum Properties {
		NUMBER_OF_PROJECTIONS;
	};

	
	// maybe it's better to move currentProjection to the viewer - in case different views can show different slices
	private int currentProjection=1;
	
	private boolean resized=false;
	
	private int width=0,height=0;
	private int originalWidth=0, originalHeight=0;
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


	private int numberOfProjections=0;
	
	private File file;

    // Actual tiltseries projections are stored inside an ImagePlus
	private ImagePlus imp=null;
	private TomoWindow window;

	// by now use a vector to store the angles
	java.util.List <Float> tiltAngles=null;
	
	// This class also saves the models of the texfields of its views, one Document per textfield
	private Document tiltTextModel=null;

	// allow Views to wait(lock) while the first (or last) image loads
	private Semaphore firstLoaded= new Semaphore(0),lastLoaded=new Semaphore(0);
	
	public TomoData(String path,TomoWindow tw){
		setFile(path);
		window=tw;
	}
	
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
	 * @return range 1..numberProjections
	 */
	public int getCurrentProjection() {
		return currentProjection;
	}

	
	/**
	 * @param currentProjection range 1..numberProjections
	 */
	public void setCurrentProjection(int currentProjection) {
		
		if((currentProjection<1) || (currentProjection > getNumberOfProjections())){
			Xmipp_Tomo.debug("setCurrentProjection("+currentProjection+")");		
			return;
		}
		
		this.currentProjection = currentProjection;
		// update all things depending on current projection/slice
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
	
	private List<Float> getTiltAngles(){
		return tiltAngles;
	}
	
	public Iterator<Float> getTiltAnglesIterator(){
		return tiltAngles.iterator();
	}
	
	/**
	 * @return null if tilt is undefined, tilt value otherwise
	 * 
	 */
	public Float getCurrentTilt(){
		Float t=null;
		try{
			t=getTiltAngle(getCurrentProjection());
		}catch (ArrayIndexOutOfBoundsException ex){
			t=null;
		}
		return t;
	}
	
	// Helper methods to manage tilt angles "list"
	
	// tiltangles list -> 0..N-1
	private Float getTiltAngle(int i){
		if(tiltAngles == null)
			return null;
		return getTiltAngles().get(i-1);
	}
	
	
	/**
	 * @param i angle index - from 1 to numberOfProjections
	 * @param t tilt angle
	 */
	private void setTiltAngle(int i, float t){
		// tiltangles list -> 0..N-1
		
		// set() below needs a non-empty vector (size > 0)
		if(tiltAngles == null){
			Vector <Float> v = new Vector<Float>();
			v.setSize(getNumberOfProjections());
			tiltAngles=v;
		}
			
		getTiltAngles().set(i-1, new Float(t));
		if(i==getCurrentProjection())
			setTiltText(t);
	}
	
	public void setTiltAnglesStep(double start, double step){
		float t=(float)start;
		for(int i=1;i<=getNumberOfProjections();i++){
			setTiltAngle(i, (float)Math.ceil(t));
			t+=step;
		}
		// setCurrentTilt(getCurrentTilt());
	}
	
	public void setTiltAngles(double start, double end){
		double step= (Math.abs(end-start) / (getNumberOfProjections() - 1));
		setTiltAnglesStep(start, step);
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
	
	public void setCurrentTilt(float t){
		setTiltAngle(getCurrentProjection(),t);
	}

	
	// method for loading the angles in sequence (it assumes an implicit order) 
	public void addTiltAngle(Float t){
		// add() below needs an empty vector
		if(tiltAngles == null)
			tiltAngles = new Vector<Float>();
			
		getTiltAngles().add(t);
		if(tiltAngles.size() == getCurrentProjection())
			setTiltText(t);
	}
	
	// prepare tilt angles array for later manual setting
	/* public void emptyTiltAngles(int count){
		tiltAngles=new Vector<Float> (count);
	}*/
	
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
	
	
	/**
	 * wrapper to simulate a settext operation using the Document methods available
	 * @param d
	 * @param t
	 */
	private void setDocumentText(Document d,String t){
		try{
			d.remove(0,d.getLength());
			d.insertString(0, t, null);
		}catch (BadLocationException e){}
	}
	
	private void setTiltText(float tilt){
		setDocumentText(getTiltModel(),String.valueOf(tilt));	
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
	
	public void loadCanceled(){
		firstImageLoaded();
	}
	
	public void waitForLastImage() throws InterruptedException{
		// Xmipp_Tomo.debug("waitForFirstImage");
		lastLoaded.acquire();
	}
	
	public void lastImageLoaded(){
		lastLoaded.release();
		// Xmipp_Tomo.debug("firstImageLoaded");
	}
	
	public void addProjection(ImageProcessor imageProcessor){
		if(getImage()==null){
			// cannot add empty stack, for example in the constructor. better init all here			
			ImageStack stack=new ImageStack(imageProcessor.getWidth(), imageProcessor.getHeight());
			stack.addSlice(null, imageProcessor);
			setImage(new ImagePlus(getFileName(), stack));
			getImage().setWindow(window);
		}else		
			getImage().getStack().addSlice(null, imageProcessor);
		
		setNumberOfProjections(getNumberOfProjections()+1);
		if(getNumberOfProjections() == 1)
			firstImageLoaded();
	}


	public boolean isResized() {
		return resized;
	}
	
	public boolean isLoadCanceled(){
		return(window.getLastCommandState() == Command.State.CANCELED);
	}


	public void setResized(boolean resized) {
		this.resized = resized;
	}
}
