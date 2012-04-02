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

import javax.swing.JPopupMenu;

import ij.ImagePlus;
import xmipp.ij.commons.ImagePlusLoader;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.viewer.ImageDimension;

public class VolumeGallery extends ImageGallery {
	protected String volFn;
	protected long volNumber;
	ImageGeneric volume;

	public VolumeGallery(GalleryData data) throws Exception {
		super(data);
		data.normalize = true; // volumes are displayed with global
								// normalization by default
		data.selection = new boolean[n];
		calculateMinAndMax();
	}

	private static final long serialVersionUID = 1L;

	// Load initial dimensions
	protected ImageDimension loadDimension() throws Exception {
		volFn = Filename.getFilename(data.selectedVolFn);
		volNumber = Filename.getNimage(data.selectedVolFn);
		volume = new ImageGeneric(data.selectedVolFn); // read image header
		volume.read(volNumber);
		if (data.resliceView != ImageGeneric.Z_NEG)
			volume.reslice(data.resliceView );
		ImageDimension dim = new ImageDimension(volume);
		return dim;
	}

//	@Override
//	protected void setZoomValue(int z) {
//		super.setZoomValue(z);
//		try {
//			volume.read(thumb_width, thumb_height, volNumber);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//	}

	@Override
	protected double[] getMinAndMax() {
		try {
			// ImageGeneric image = new ImageGeneric(volFn);
			// image.read(volNumber);
			double[] stats = volume.getStatistics();
			// image.destroy();
			return stats;
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return null;
	}

	@Override
	public String getItemKey(int index) throws Exception {
		return String.format("%d_(%d,%d)", index, thumb_width, thumb_height);
	}

	@Override
	public String getTitle() {
		return String.format("Volume: %s (%d x %d x %d)", data.selectedVolFn,
				image_width, image_height, n);
	}

	@Override
	protected ImageItem createItem(int index, String key) throws Exception {
		ImageGeneric preview = new ImageGeneric();
		volume.getPreview(preview, thumb_width, thumb_height, 
				index + 1, ImageGeneric.FIRST_IMAGE);
		ImagePlus imp = XmippImageConverter.convertToImagePlus(preview);
		ImageItem item = new ImageItem(index);
		item.setImagePlus(imp);
		return item;
	}

	@Override
	public boolean handleDoubleClick(int row, int col) {
		try {
			int index = getIndex(row, col);
			ImagePlus imp = XmippImageConverter.convertToImagePlus(volume,
					ImageGeneric.FIRST_IMAGE, index + 1);
			new XmippImageWindow(new ImagePlusLoader(imp), getLabel(row, col));
			return true;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return false;
	}// function handleDoubleClick

	@Override
	public String getLabel(int row, int col) {
		int index = getIndex(row, col);
		return String.format("slice %d", index + 1);
	}

	@Override
	public ImagePlus getImagePlus() {
		try {
			return XmippImageConverter.convertToImagePlus(volume);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	@Override
	public boolean handleRightClick(int row, int col,
			XmippPopupMenuCreator xpopup) {
		return false;
	}

	// @Override
	// public String getImageFilenameAt(int row, int col) {
	// return null;
	// }

}// class VolumeGallery
