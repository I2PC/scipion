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

import ij.ImagePlus;

import java.io.File;
import java.util.ArrayList;

import xmipp.ij.commons.ImagePlusLoader;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.utils.DEBUG;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.viewer.ImageDimension;

public class MetadataGallery extends ImageGallery {
	private static final long serialVersionUID = 1L;

	// Label to be rendered
	protected ColumnInfo renderLabel;
	protected ColumnInfo displayLabel;
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
		// data.globalRender = false;

		for (ColumnInfo ci : data.labels) {
			if (ci.visible)
				visibleLabels.add(ci);
			// if (ci.render)
			// data.globalRender = true;
		}
		ImageDimension dim = null;

		if (data.hasRenderLabel()) {
			renderLabel = data.ciFirstRender;
			displayLabel = renderLabel;
			// if (renderLabels) {
			for (int i = 0; i < data.ids.length; ++i) {
				String imageFn = getImageFilename(i, renderLabel.getLabel());
				if (imageFn != null && Filename.exists(imageFn)) {
					try {
					image = new ImageGeneric(imageFn);
					dim = new ImageDimension(image);
					break;
					} catch (Exception e){
						dim = null;
					}
				}
			}
		}

		if (dim == null)
			dim = new ImageDimension(30);
		dim.setZDim(data.ids.length);

		return dim;
	}

	@Override
	protected ImageItem createItem(int index, String key) throws Exception {
		return createImageItem(index, renderLabel.getLabel(),
				displayLabel.getLabel(), key);
	}

	public String getLabel(int row, int col) {
		try {
			int index = getIndex(row, col);
			long objId = data.ids[index];
			if (data.is2dClassification){
				int ref = data.md.getValueInt(MDLabel.MDL_REF, objId);
				long count = data.md.getValueLong(MDLabel.MDL_CLASS_COUNT, objId);
				return String.format("class %d (%d images)", ref, count);
			}
			else
				return data.md.getValueString(displayLabel.getLabel(), objId);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * Function to create an image item
	 */
	protected ImageItem createImageItem(int index, int renderLabel,
			int displayLabel, String key) throws Exception {
		String imageFn = getImageFilename(index, renderLabel);
		long objId = data.ids[index];
		ImageItem item = new ImageItem(index);
		ImagePlus imp = null;
		if (imageFn != null && Filename.exists(imageFn)) {
			if (data.useGeo)
				image.readApplyGeo(imageFn, data.md, objId, thumb_width,
						thumb_height, data.wrap);
			else
				image.read(imageFn, thumb_width, thumb_height);
			imp = XmippImageConverter.convertToImagePlus(image);
		}
		item.setImagePlus(imp);
		return item;
	}

	@Override
	public String getItemKey(int index) throws Exception {
		return getItemKey(index, renderLabel.getLabel());
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
		// String key = String.format(format, thumb_width, thumb_height);
		// DEBUG.printMessage(String.format("key: %s", key));
		return String.format(format, thumb_width, thumb_height);
	}

	@Override
	public String getTitle() {
		String title = "Metadata: " + (filename != null ? 
				Filename.getBaseName(filename) : "");
		if (n > 1)
			title += String.format(" %d items", n);
		if (data.hasRenderLabel())
			title += String.format(" (%d x %d)", image_width, image_height);
		return title;
	}//function getTitle

	public String getImageFilename(int index, int label) {
		try {
			return data.getValueFromLabel(index, label);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}//function getImageFilename

//	@Override
//	public String getImageFilenameAt(int row, int col) {
//		return data.isImageFile(renderLabel) ? data.getValueFromCol(row,
//				renderLabel) : null;
//	}
	
	@Override
	public boolean handleDoubleClick(int row, int col){
		try {
			if (data.isImageFile(renderLabel)) {
				new XmippImageWindow(new ImagePlusLoader(data.getValueFromCol(getIndex(row, col), renderLabel)));
				return true;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return false;
	}//function handleDoubleClick

	@Override
	protected double[] getMinAndMax() {
		try {
			return data.md.getStatistics(false);
		} catch (Exception ex) {
			DEBUG.printException(ex);
		}
		return null;
	}//function getMinAndMax

	/** Change the use of geometry info */
	public void setUseGeometry(boolean geo, boolean wrap) {
		if (!geo)
			wrap = false;
		boolean changed = data.useGeo != geo || data.wrap != wrap;
		data.useGeo = geo;
		data.wrap = wrap;
		if (changed)
			fireTableDataChanged();
	}//function setUseGeometry

	@Override
	public ImagePlus getImagePlus() {
		try {
			return XmippImageConverter.readMetadataToImagePlus(data.md);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}//function getImagePlus

	@Override
	public boolean handleRightClick(int row, int col,
			XmippPopupMenuCreator xpopup) {
		return true;
	}//function handleRightClick
}
