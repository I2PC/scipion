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

import java.util.ArrayList;

import ij.ImagePlus;
import xmipp.ij.XmippImageConverter;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;

public class MetadataGallery extends ImageGallery {
	private static final long serialVersionUID = 1L;

	// Label to be rendered
	protected int renderLabel;
	protected int displayLabel;

	// Also store the visible ones to fast access
	ArrayList<ColumnInfo> visibleLabels;

	public MetadataGallery(GalleryData data) throws Exception {
		super(data);
		data.normalize = false;
	}

	/** Update the columns display information */
	public void updateColumnInfo(ArrayList<ColumnInfo> newInfo) {
		int n = newInfo.size();
		boolean changed = false;
		data.globalRender = false;

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j) {
				ColumnInfo ci1 = data.labels.get(i);
				ColumnInfo ci2 = newInfo.get(j);
				if (ci1.label == ci2.label) {
					if (ci1.updateInfo(ci2))
						changed = true;

					if (ci1.visible)
						visibleLabels.add(ci1);
					if (ci1.render)
						data.globalRender = true;
					if (i != j)
						changed = true;
				}
			}

		if (changed) {
			data.labels = newInfo;
			visibleLabels.clear();
			for (ColumnInfo ci : data.labels)
				if (ci.visible)
					visibleLabels.add(ci);
			calculateCellSize();
			cols = visibleLabels.size();
			// fireTableDataChanged();
			fireTableStructureChanged();
		}
	}

	// Load initial dimensions
	protected ImageDimension loadDimension() throws Exception {
		// Set information about columns
		visibleLabels = new ArrayList<ColumnInfo>();
		data.globalRender = false;

		for (ColumnInfo ci : data.labels) {
			if (ci.visible)
				visibleLabels.add(ci);
			if (ci.render)
				data.globalRender = true;
		}
		ImageDimension dim = null;

		if (data.hasRenderLabel()) {
			renderLabel = data.ciFirstRender.label;
			displayLabel = renderLabel;
			// if (renderLabels) {
			for (int i = 0; i < data.ids.length; ++i) {
				ImageGeneric image = getImage(i, renderLabel);
				if (image != null) {
					dim = new ImageDimension(image);
					image.destroy();
					break;
				}
			}
		}

		if (dim == null)
			dim = new ImageDimension(30, 30, 0, 0);
		dim.setZDim(data.ids.length);

		return dim;
	}

	@Override
	protected ImageItem createItem(int index, String key) throws Exception {
		return createImageItem(index, renderLabel, displayLabel, key);
	}

	/**
	 * Function to create an image item
	 * 
	 * @throws Exception
	 */
	protected ImageItem createImageItem(int index, int renderLabel,
			int displayLabel, String key) throws Exception {
		ImageGeneric image = getImage(index, renderLabel);
		String labelStr = data.md.getValueString(displayLabel, data.ids[index]);
		ImagePlus imp = null;
		if (image != null) {
			if (data.useGeo)
				image.readApplyGeo(data.md, data.ids[index], thumb_width,
						thumb_height);
			else
				image.read(thumb_width, thumb_height);
			imp = XmippImageConverter.convertImageGenericToImageJ(image);
			
		}
		return new ImageItem(key, labelStr, imp);
	}

	@Override
	protected String getItemKey(int index) throws Exception {
		return getItemKey(index, renderLabel);
	}

	/**
	 * Return a key string using label
	 * 
	 * @throws Exception
	 */
	protected String getItemKey(int index, int label) throws Exception {
		String format = data.md.getValueString(label, data.ids[index])
				+ "(%d,%d)";
		if (data.useGeo)
			format += "_geo";
		String key = String.format(format, thumb_width, thumb_height);
		// DEBUG.printMessage(String.format("key: %s", key));
		return String.format(format, thumb_width, thumb_height);
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
	protected ImageGeneric getImage(int index, int label) {
		try {
			String imgFn = data.md.getValueString(label, data.ids[index]);
			if (Filename.exists(imgFn)){
				ImageGeneric image = new ImageGeneric(imgFn);
				return image;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	@Override
	protected double[] getMinAndMax() {
		try {
			return data.md.getStatistics(false);
		} catch (Exception ex) {
			DEBUG.printException(ex);
		}
		return null;
	}

	/** Change the use of geometry info */
	public void setUseGeometry(boolean value) {
		if (data.useGeo != value) {
			data.useGeo = value;
			fireTableDataChanged();
		}
	}

	@Override
	public ImagePlus getImagePlus() {
		try {
			return XmippImageConverter.readMetadataToImageJ(data.md);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
}
