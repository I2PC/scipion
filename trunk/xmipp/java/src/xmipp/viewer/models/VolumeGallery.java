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
import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.viewer.ImageDimension;
import xmipp.viewer.ImageItem;

public class VolumeGallery extends ImageGallery {
	protected String volFn;
	protected long volNumber;
	ImageGeneric volHeader;
	
	public VolumeGallery(GalleryData data) throws Exception {
		super(data);
		data.normalize = true; // volumes are displayed with global normalization by default
		volFn = Filename.getFilename(data.selectedVol);
		volNumber = Filename.getNimage(data.selectedVol);
		data.selection = new boolean[n];
		calculateMinAndMax();
	}

	private static final long serialVersionUID = 1L;

	// Load initial dimensions
	protected ImageDimension loadDimension() throws Exception {
		volHeader = new ImageGeneric(data.selectedVol); // read image header
		ImageDimension dim = new ImageDimension(volHeader);
		return dim;
	}
	
	@Override
	protected void setZoomValue(int z) {
		super.setZoomValue(z);
		try {
			volHeader.read(thumb_width, thumb_height, volNumber);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	@Override
	protected double[] getMinAndMax() {
		try {
			ImageGeneric image = new ImageGeneric(volFn);
			image.read(volNumber);
			double[] stats = image.getStatistics();
			image.destroy();
			return stats;
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return null;
	}

	@Override
	protected String getItemKey(int index) throws Exception {
		return String.format("%d_(%d,%d)", index, thumb_width, thumb_height);
	}
	
	@Override
	public String getTitle() {
		return String.format("Volume: %s (%d x %d x %d)", data.selectedVol, image_width, image_height, n);
	}

	@Override
	protected ImageItem createItem(int index, String key) throws Exception {
		ImagePlus imp = XmippImageConverter.convertToImagePlus(volHeader, ImageGeneric.FIRST_IMAGE, index + 1);
		String label = String.format("%d", index + 1);
		return new ImageItem(key, label, imp);
	}

	@Override
	public ImagePlus getImagePlus() {
		try {
			return XmippImageConverter.readToImagePlus(volHeader);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	@Override
	public boolean handleRightClick(int row, int col, XmippPopupMenuCreator xpopup) {
		return false;
	}

}
