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
import java.io.FileNotFoundException;
import java.io.IOException;

import ij.ImagePlus;
import xmipp.ImageDouble;
import xmipp.MetaData;

/**
 * Why? Model the active stack that will be displayed, isolating the viewer from the
 * gory details (current slice, current stack, cache, etc)
 */
public class StackModel extends AbstractModel{
	// each file is accessed through its MetaData, hence somewhere every MetaData in use
	// must be kept as the workflow progresses. For example, keep in the Workflow model
	// the filepath and MetaData of each node
	
	private Cache <String,ImageDouble> imageDoubleCache;
	private Cache <String,ImagePlus> imagePlusCache;
	
	private Workflow workflow;
	
	public static enum Properties {
		NUMBER_OF_PROJECTIONS,CURRENT_PROJECTION_NUMBER,CURRENT_TILT_ANGLE,CURRENT_PROJECTION_ENABLED, CURRENT_PROJECTION;
	};
	
	// maybe it's better to move currentProjection to the viewer - in case different views can show different slices
	// TODO: strictly speaking, the order of the projections is not specified by the sequence of the images in the
	// stack file, it is specified by the tilt angle. It's just a coincidence that the first image in the sequence is
	// the one with the smallest tilt angle.
	private int currentProjectionNumber=1;
	private int numberOfProjections=0;
	
	public StackModel(Workflow workflow){
		this.workflow=workflow;
	}
	
	private Workflow getWorkflow() {
		return workflow;
	}


	private Cache <String,ImageDouble> getImageDoubleCache() {
		if(imageDoubleCache == null)
			imageDoubleCache = new Cache<String,ImageDouble>();
		return imageDoubleCache;
	}

	private Cache <String,ImagePlus> getImagePlusCache() {
		if(imagePlusCache == null)
			imagePlusCache = new Cache<String,ImagePlus>();
		return imagePlusCache;
	}
	
	private String getCurrentImageFilePath(){
		String ret = null;
		UserActionIO currentIo=getCurrentIO();
		if(currentIo==null)
			return ret;
		
		ret = currentIo.getFilePath(getCurrentProjectionNumber());
		
		return ret;
	}
	
	public ImagePlus getCurrentImage(){
		ImagePlus ret = new ImagePlus();

		String filePath = getCurrentImageFilePath();
		if(filePath == null)
			return ret;
		
		ret= getImagePlusCache().get(filePath);
		if(ret != null)
			return ret;
		
		ImageDouble image = getImageDoubleCache().get(filePath);
		if(image == null){
			// cache miss
			image = new ImageDouble();
			try{
				image.readSlice(filePath);
			} catch (FileNotFoundException ex) {
				Logger.debug("loadEMBackground - ", ex);
			} catch (IOException ex) {
				Logger.debug("loadEMBackground - Error opening file ",
						ex);
			} catch (OutOfMemoryError err) {
				Logger.debug("loadEMBackground - Out of memory"
						+ err.toString());
			} catch (Exception ex) {
				Logger.debug("loadEMBackground - unexpected exception",
						ex);
			}
			getImageDoubleCache().put(filePath, image);
		}
		ret= Converter.convertToImagePlus(image);
		getImagePlusCache().put(filePath, ret);
		return ret;
	}
	
	public void updateCurrentImage(ImagePlus image){
		updateImage(getCurrentImageFilePath(),image);
	}
	
	public void updateImage(String filePath,ImagePlus image){
		if( (filePath == null) || (image == null))
			return;
		if(getImagePlusCache().get(filePath) == null)
			return;
		getImagePlusCache().put(filePath, image);
		if(filePath.equals(getCurrentImageFilePath()))
			firePropertyChange(Properties.CURRENT_PROJECTION.name(), false, true);
	}

	public void updateImage(String filePath,ImageDouble image){
		updateImage(filePath,image,null);
	}
	
	public void updateImage(String filePath,ImageDouble image, ImagePlus ip){
		if( (filePath == null) || (image == null))
			return;
		if(getImageDoubleCache().get(filePath) != null)
			getImageDoubleCache().put(filePath, image);
		
		if (getImagePlusCache().get(filePath) != null){
			if(ip == null)
				ip= Converter.convertToImagePlus(image);
			updateImage(filePath,ip);
		}
	}
	
	private UserActionIO getCurrentIO(){
		return getWorkflow().getCurrentUserActionIO();
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
			Logger.debug("setCurrentProjection("+newProjectionNumber+") - out of bounds");		
			return;
		}
		int oldProjectionNumber=currentProjectionNumber;
		this.currentProjectionNumber = newProjectionNumber;
	// possible dependencies on currentProjectionNumber should be handled by the View, not by this model
		// Unluckily, ImagePlus mixes data and presentation, so it's better to keep ImagePlus related operations here
		// TODO: update image displayed?
		/*
		try{
			getImage().setSlice(getCurrentProjectionNumber());
		} catch (IllegalArgumentException ex){
			Logger.debug("setCurrentProjection - ", ex);
		}*/
		firePropertyChange(Properties.CURRENT_PROJECTION_NUMBER.name(), oldProjectionNumber, currentProjectionNumber);
	}
	
	/**
	 * @return range 0..N
	 */
	public int getNumberOfProjections(){
		// get the number from the metadata
		return getCurrentIO().getNumberOfProjections();
	}
	
	private void setNumberOfProjections(int n){
		int oldNumberOfProjections=numberOfProjections;
		numberOfProjections=n;
		if(n > 1)
			firePropertyChange(Properties.NUMBER_OF_PROJECTIONS.name(), oldNumberOfProjections, n);

	}
	
	// TODO: complete (from Tomodata)
	public boolean isCurrentEnabled(){
		return true;
	}
	
	/**
	 * @return 1..N, cyclic (next to last is 1)
	 */
	public int getNextProjection(){
		int ret= getCurrentProjectionNumber();
		try{
			ret= (getCurrentProjectionNumber() % getNumberOfProjections()) + 1;
		}catch (ArithmeticException ex){
			Logger.debug("Error: ", ex);
		}
		return ret;
	}
	
}
