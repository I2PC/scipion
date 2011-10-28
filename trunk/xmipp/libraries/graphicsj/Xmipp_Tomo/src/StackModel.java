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
import java.awt.Rectangle;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Enumeration;
import java.util.Vector;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import xmipp.ImageDouble;
import xmipp.MetaData;
import xmipp.MDLabel;

/**
 * Why? Model the active stack that will be displayed, isolating the viewer from
 * the gory details (current slice, current stack, cache, etc)
 */
public class StackModel extends AbstractModel {
	// each file is accessed through its MetaData, hence somewhere every
	// MetaData in use
	// must be kept as the workflow progresses. For example, keep in the
	// Workflow model
	// the filepath and MetaData of each node

	private Cache<String, ImageDouble> imageDoubleCache;
	private Cache<String, ImagePlus> imagePlusCache;

	private Workflow workflow;

	private double min = Double.POSITIVE_INFINITY,
			max = Double.NEGATIVE_INFINITY;

	public static enum Properties {
		NUMBER_OF_PROJECTIONS, CURRENT_PROJECTION_NUMBER, CURRENT_TILT_ANGLE, CURRENT_PROJECTION_ENABLED, CURRENT_PROJECTION;
	};

	// maybe it's better to move currentProjection to the viewer - in case
	// different views can show different slices
	// TODO: strictly speaking, the order of the projections is not specified by
	// the sequence of the images in the
	// stack file, it is specified by the tilt angle. It's just a coincidence
	// that the first image in the sequence is
	// the one with the smallest tilt angle.
	private int currentProjectionNumber = 1;
	private int numberOfProjections = 0;

	public StackModel(Workflow workflow) {
		this.workflow = workflow;
	}

	private Workflow getWorkflow() {
		return workflow;
	}

	private Cache<String, ImageDouble> getImageDoubleCache() {
		if (imageDoubleCache == null)
			imageDoubleCache = new Cache<String, ImageDouble>();
		return imageDoubleCache;
	}

	private Cache<String, ImagePlus> getImagePlusCache() {
		if (imagePlusCache == null)
			imagePlusCache = new Cache<String, ImagePlus>();
		return imagePlusCache;
	}

	private String getCurrentImageFilePath() {
		String ret = null;
		UserAction currentAction = getCurrentUserAction();
		if (currentAction == null)
			return ret;

		ret = currentAction.getOutputFilePath(getCurrentProjectionNumber());

		return ret;
	}

	private ImageDouble getCurrentImageDouble() {
		String filePath = getCurrentImageFilePath();
		if (filePath == null)
			return null;
		return getImageDoubleCache().get(filePath);
	}

	public ImagePlus getCurrentImage() {
		ImagePlus ret = new ImagePlus();

		String filePath = getCurrentImageFilePath();
		// Logger.debug(filePath);
		if (filePath == null)
			return ret;

		ret = getImagePlusCache().get(filePath);
		if (ret != null)
			return ret;

		ImageDouble image = getCurrentImageDouble();
		if (image == null) {
			// cache miss
			image = new ImageDouble();
			try {
				image.readSlice(filePath);
			} catch (FileNotFoundException ex) {
				Logger.debug("File not found ", ex);
				return null;
			} catch (IOException ex) {
				Logger.debug("Error opening file ", ex);
				return null;
			} catch (OutOfMemoryError err) {
				Logger.debug("Out of memory" + err.toString());
				return null;
			} catch (Exception ex) {
				Logger.debug("unexpected exception", ex);
				return null;
			}
			getImageDoubleCache().put(filePath, image);
		}
		ret = Converter.convertToImagePlus(image);
		getImagePlusCache().put(filePath, ret);
		// WindowManager.setTempCurrentImage(ret);
		// scale image if needed
/*		if(ret.getWidth() > 512 || ret.getHeight() > 512){
			// TODO: take into consideration height too (for non-square images)
			double mag=(512.0/ret.getWidth()) * 100;
			//IJ.run("Set... ", "zoom=" + ((int) mag));
		}*/
		return ret;
	}

	public void refreshCurrentImage() {
		// invalidate current image in ImagePlus cache, so it will be refreshed
		// on next request
		getImagePlusCache().remove(getCurrentImageFilePath());
	}

	public void updateCurrentImage(ImagePlus image) {
		updateImage(getCurrentImageFilePath(), image);
	}

	public void updateImage(String filePath, ImagePlus image) {
		if ((filePath == null) || (image == null))
			return;
		if (getImagePlusCache().get(filePath) == null)
			return;
		getImagePlusCache().put(filePath, image);
		if (filePath.equals(getCurrentImageFilePath()))
			firePropertyChange(Properties.CURRENT_PROJECTION.name(), false,
					true);
	}

	public void updateImage(String filePath, ImageDouble image) {
		updateImage(filePath, image, null);
	}

	public void updateImage(String filePath, ImageDouble image, ImagePlus ip) {
		if ((filePath == null) || (image == null))
			return;
		if (getImageDoubleCache().get(filePath) != null)
			getImageDoubleCache().put(filePath, image);

		if (getImagePlusCache().get(filePath) != null) {
			if (ip == null)
				ip = Converter.convertToImagePlus(image);
			updateImage(filePath, ip);
		}
	}

	private UserAction getCurrentUserAction() {
		return getWorkflow().getCurrentUserAction();
	}

	public String getCurrentFileName() {
		return getCurrentUserAction().getOutputFileName();
	}

	public int getCurrentWidth() {
		if (getCurrentImageDouble() == null)
			return -1;
		return getCurrentImageDouble().getXsize();
	}

	public int getCurrentHeight() {
		if (getCurrentImageDouble() == null)
			return -1;
		return getCurrentImageDouble().getYsize();
	}

	public int getCurrentBitDepth() {
		// TODO: maybe it's better getting the bit depth from the ImageDouble
		if (getCurrentImage() == null)
			return -1;
		return getCurrentImage().getBitDepth();
	}

	/**
	 * @return range 1..numberProjections
	 */
	public int getCurrentProjectionNumber() {
		return currentProjectionNumber;
	}

	public String getCurrentProjectionInfo() {
		return "Projection #" + getCurrentProjectionNumber() + " ("
				+ getCurrentUserActionInfo();
	}

	private String getCurrentUserActionInfo() {
		return getCurrentUserAction().getInfo(getCurrentProjectionNumber());
	}

	private long getCurrentProjectionId() {
		return getCurrentUserAction().getProjectionId(
				getCurrentProjectionNumber());
	}

	/**
	 * @param newProjectionNumber
	 *            range 1..numberProjections
	 */
	public void setCurrentProjectionNumber(int newProjectionNumber) {

		if ((newProjectionNumber < 1)
				|| (newProjectionNumber > getNumberOfProjections())) {
			Logger.debug("setCurrentProjection(" + newProjectionNumber
					+ ") - out of bounds");
			return;
		}
		int oldProjectionNumber = currentProjectionNumber;
		this.currentProjectionNumber = newProjectionNumber;

		firePropertyChange(Properties.CURRENT_PROJECTION_NUMBER.name(),
				oldProjectionNumber, currentProjectionNumber);
	}

	/**
	 * @return range 0..N
	 */
	public int getNumberOfProjections() {
		// get the number from the metadata
		return getCurrentUserAction().getNumberOfProjections();
	}

	private void setNumberOfProjections(int n) {
		int oldNumberOfProjections = numberOfProjections;
		numberOfProjections = n;
		if (n > 1)
			firePropertyChange(Properties.NUMBER_OF_PROJECTIONS.name(),
					oldNumberOfProjections, n);

	}

	/**
	 * @return 1..N, cyclic (next to last is 1)
	 */
	public int getNextProjection() {
		int ret = getCurrentProjectionNumber();
		try {
			ret = (getCurrentProjectionNumber() % getNumberOfProjections()) + 1;
		} catch (ArithmeticException ex) {
			Logger.debug("Error: ", ex);
		}
		return ret;
	}

	public void discardCurrentProjection() {
		if (getNumberOfProjections() > 1) {
			setEnabled(getCurrentProjectionId(), false);

		}

	}
	
	public void enableCurrentProjection() {
		if (getNumberOfProjections() > 1) {
			setEnabled(getCurrentProjectionId(), true);
		}
	}

	private void setEnabled(long id, boolean e) {
		// getEnabledProjections().set(i-1, new Boolean(e));
		int enabled = 0;
		if (e)
			enabled = 1;

		getCurrentUserAction().setEnabled(id, enabled);
		firePropertyChange(Properties.CURRENT_PROJECTION_ENABLED.name(), true,false);
	}

	// Tilt features

	public Double getCurrentTiltAngle() {
		Double t = null;
		try {
			t = getTiltAngle(getCurrentProjectionNumber());
		} catch (ArrayIndexOutOfBoundsException ex) {
			t = null;
		}
		return t;
	}

	public Enumeration<Double> getTiltAngles() {
		int size = getNumberOfProjections();

		Vector<Double> tiltAngles = new Vector<Double>(size);
		int i;
		for (i = 1; i < size + 1; i++)
			tiltAngles.add(getTiltAngle(i));
		return tiltAngles.elements();
	}

	/**
	 * @param i
	 *            1..N
	 */
	private Double getTiltAngle(long id) {
		return getCurrentUserAction().getTiltAngle(id);
	}

	// right now return only grayscale value as double
	// more complete method at ImagePlus.getValueAsString()
	public double getPixelValue(int x, int y) {
		int[] v = getCurrentImage().getPixel(x, y);
		// 32bit images
		return Float.intBitsToFloat(v[0]);
	}

	public boolean isCurrentEnabled() {
		return getCurrentUserAction().isEnabled(getCurrentProjectionId());
	}


	// TODO: Min an Max not set. We need to load all the projections of the stack to know the min and max
	public void normalize() {
		getCurrentImage().getProcessor().setMinAndMax(getMin(), getMax());
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

	void applySelFile() {
		getCurrentUserAction().applySelFile();
		
		/*
		 * String outputFile = "prueba.sel"; int size=getNumberOfProjections();
		 * FileWriter fw = null; PrintWriter pw = null; int enable; try { fw =
		 * new FileWriter(outputFile); pw = new PrintWriter(fw);
		 * 
		 * pw.println(
		 * "# XMIPP_STAR_1 *\n#\ndata_\nloop_\n_image\n_enable\n_angleTilt");
		 * for (int i = 1; i < size+1; i++){ if (isEnabled(i)) enable=1; else
		 * enable=0; pw.println(getImageFilePath(i) + "\t\t" +enable +
		 * "\t"+getTiltAngle(i) ); } } catch (Exception e) {
		 * e.printStackTrace(); } finally { try {
		 * 
		 * if (null != fw) fw.close(); } catch (Exception e2) {
		 * 
		 * }
		 */ 
	}

}
