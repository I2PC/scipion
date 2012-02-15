/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
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

package xmipp.viewer;

import ij.ImagePlus;
import xmipp.ij.XmippImageConverter;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;

public class MetadataGallery extends ImageGallery {

	private static final long serialVersionUID = 1L;

	// MetaData with elements to display
	protected MetaData md;
	// Label to be rendered
	protected int renderLabel;
	protected int displayLabel;
	
	// Store labels and names
	int[] labels;
	String[] labelStr;

	// Ids of objects stored in metadata
	long[] ids;

	public MetadataGallery(String fn, int zoom) throws Exception {
		super(fn, zoom);
		labels = md.getActiveLabels();
		labelStr = new String[labels.length];
		for (int i = 0; i < labels.length; ++i)
			labelStr[i] = MetaData.label2Str(labels[i]);
	}

	// Load initial dimensions
	protected ImageDimension loadDimension() throws Exception {
		md = new MetaData(filename);
		ids = md.findObjects();
		renderLabel = MDLabel.MDL_IMAGE;
		displayLabel = MDLabel.MDL_IMAGE;
		ImageGeneric image = getImage(0, renderLabel);
		ImageDimension dim = new ImageDimension(image);
		// TODO: check this well, now asuming not volumes in metadata
		dim.setZDim(ids.length);
		image.destroy();
		return dim;
	}
	
	@Override
	protected ImageItem createItem(int index, String key) throws Exception {
		return createImageItem(index, renderLabel, displayLabel, key);
	}
	
	/** Function to create an image item 
	 * @throws Exception */
	protected ImageItem createImageItem(int index, int renderLabel, int displayLabel, String key) throws Exception{
		ImageGeneric image = getImage(index, renderLabel);
		image.readApplyGeo(md, ids[index], thumb_width, thumb_height);
		ImagePlus imp = XmippImageConverter.convertImageGenericToImageJ(image);
		String labelStr = md.getValueString(displayLabel, ids[index]);
		return new ImageItem(key, labelStr, imp);		
	}

	@Override
	protected String getItemKey(int index) throws Exception {
		return getItemKey(index, renderLabel);
	}
	
	/** Return a key string using label 
	 * @throws Exception */
	protected String getItemKey(int index, int label) throws Exception{
		String filename = md.getValueString(label, ids[index]);
		return String.format("%s(%d,%d)", filename, thumb_width, thumb_height);
	}
	
	@Override
	public String getTitle() {
		return String.format("Metadata: %s (%d)", filename, n);
	}

	/**
	 * Create and returns the image generic from a given an index
	 * 
	 * @param index
	 *            Index of the image
	 * @return ImageGeneric created
	 * @throws Exception
	 *             if can not load image
	 */
	protected ImageGeneric getImage(int index, int label) throws Exception {
		String imgFn = md.getValueString(label, ids[index]);
		return new ImageGeneric(imgFn);
	}

	@Override
	protected double[] getMinAndMax() {
		try {
			return md.getStatistics(false);
		} catch (Exception ex) {
			DEBUG.printException(ex);
		}
		return null;
	}
}
