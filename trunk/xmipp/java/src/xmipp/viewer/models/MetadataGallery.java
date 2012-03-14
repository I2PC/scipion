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

package xmipp.viewer.models;

import java.util.ArrayList;

import javax.swing.JPopupMenu;

import ij.ImagePlus;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.viewer.ImageDimension;
import xmipp.viewer.ImageItem;

public class MetadataGallery extends ImageGallery {
	private static final long serialVersionUID = 1L;

	// Label to be rendered
	protected int renderLabel;
	protected int displayLabel;
	protected ImageGeneric image;

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
	@Override
	protected ImageDimension loadDimension() throws Exception {
		// Set information about columns
		visibleLabels = new ArrayList<ColumnInfo>();
//		data.globalRender = false;

		for (ColumnInfo ci : data.labels) {
			if (ci.visible)
				visibleLabels.add(ci);
//			if (ci.render)
//				data.globalRender = true;
		}
		ImageDimension dim = null;

		if (data.hasRenderLabel()) {
			renderLabel = data.ciFirstRender.label;
			displayLabel = renderLabel;
			// if (renderLabels) {
			for (int i = 0; i < data.ids.length; ++i) {
				String imageFn = getImageFilename(i, renderLabel);
				if (imageFn != null && Filename.exists(imageFn)) {
					image = new ImageGeneric(imageFn); 
					dim = new ImageDimension(image);
					//image.destroy();
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
	 */
	protected ImageItem createImageItem(int index, int renderLabel,
			int displayLabel, String key) throws Exception {
		String imageFn = getImageFilename(index, renderLabel);
		long objId = data.ids[index];
		String labelStr = data.md.getValueString(displayLabel, objId);
		ImagePlus imp = null;
		if (imageFn != null) {
			if (data.useGeo)
				image.readApplyGeo(imageFn, data.md, objId, thumb_width,
						thumb_height, data.wrap);
			else
				image.read(imageFn, thumb_width, thumb_height);
			imp = XmippImageConverter.convertToImagePlus(image);
			//image.destroy();
			
		}
		ImageItem ii = new ImageItem(key, labelStr, imp);
		ii.isEnabled = data.md.getEnabled(objId);
		return ii;
	}

	@Override
	protected String getItemKey(int index) throws Exception {
		return getItemKey(index, renderLabel);
	}

	/**
	 * Return a key string using label
	 */
	protected String getItemKey(int index, int label) throws Exception {
		String format = data.md.getValueString(label, data.ids[index])
				+ "(%d,%d)";
		if (data.useGeo)
			format += "_geo";
		if (data.wrap)
			format += "_wrap";
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
	public String getImageFilename(int index, int label) {
		try {
			return data.getValueFromLabel(index, label);
//			if (Filename.exists(imgFn)){
//				ImageGeneric image = new ImageGeneric(imgFn);
//				return image;
//			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	@Override
	public String getImageFilenameAt(int row, int col){
		if (data.isImageFile(col))
			return data.getValueFromCol(row, col);
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
	public void setUseGeometry(boolean geo, boolean wrap) {
		if (!geo)
			wrap = false;
		boolean changed = data.useGeo != geo || data.wrap != wrap;
		data.useGeo = geo;
		data.wrap = wrap;
		if (changed) 
			fireTableDataChanged();
	}

	@Override
	public ImagePlus getImagePlus() {
		try {
			return XmippImageConverter.readMetadataToImagePlus(data.md);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	@Override
	public boolean handleRightClick(int row, int col, XmippPopupMenuCreator xpopup) {
		return true;
	}
}
