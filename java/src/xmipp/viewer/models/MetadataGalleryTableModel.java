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
import xmipp.viewer.windows.ImagesWindowFactory;

public class MetadataGalleryTableModel extends ImageGalleryTableModel
{

	private static final long serialVersionUID = 1L;

	// Label to be rendered
	protected ColumnInfo renderLabel;
	protected ColumnInfo displayLabel;
	protected ImageGeneric image;

	// Also store the visible ones to fast access
	public ArrayList<ColumnInfo> visibleLabels;

	public MetadataGalleryTableModel(GalleryData data) throws Exception
	{
		super(data);
		data.normalize = false;

	}

	/** Update the columns display information */
	public void updateColumnInfo(ArrayList<ColumnInfo> newInfo)
	{
		int n = newInfo.size();
		boolean changed = false;
		data.globalRender = false;
		data.ciFirstRender = null;

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
			{
				ColumnInfo ci1 = data.labels.get(i);
				ColumnInfo ci2 = newInfo.get(j);
				if (ci1.label == ci2.label)
				{
					if (ci1.updateInfo(ci2))
						changed = true;

					if (ci1.visible)
					{
						visibleLabels.add(ci1);
						if (data.ciFirstRender == null && ci1.allowRender)
							data.ciFirstRender = ci1;
					}
					if (ci1.render)
						data.globalRender = true;
					if (i != j)
						changed = true;
				}
			}

		if (changed)
		{
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
	
	

	public ImagePlus getImage(long id, String imagepath, int width, int height, boolean useGeo, boolean wrap)
	{
		ImagePlus imp = null;
		if (imagepath != null && Filename.exists(imagepath))
		{
			try
			{
				imp = XmippImageConverter.readMdRowToImagePlus(imagepath, data.md, id, width, height, useGeo, wrap);
			}
			catch (Exception ex)
			{
				imp = null;
			}
		}
		return imp;
	}

	// Load initial dimensions
	@Override
	protected ImageDimension loadDimension() throws Exception
	{
		// Set information about columns
		visibleLabels = new ArrayList<ColumnInfo>();
		// data.globalRender = false;

		for (ColumnInfo ci : data.labels)
		{
			if (ci.visible)
				visibleLabels.add(ci);
			// if (ci.render)
			// data.globalRender = true;
		}
		ImageDimension dim = null;

		if (data.hasRenderLabel())
		{
			renderLabel = data.ciFirstRender;
			displayLabel = renderLabel;
			// if (renderLabels) {
			for (int i = 0; i < data.ids.length; ++i)
			{
				String imageFn = getImageFilename(i, renderLabel.getLabel());
				if (imageFn != null && Filename.exists(imageFn))
				{
					try
					{
						image = new ImageGeneric(imageFn);
						dim = new ImageDimension(image);
						break;
					}
					catch (Exception e)
					{
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
	protected ImageItem createItem(int index, String key) throws Exception
	{
		return createImageItem(index, renderLabel.getLabel(), displayLabel.getLabel(), key);
	}

	public String getLabel(int row, int col)
	{
		try
		{
			int index = getIndex(row, col);
			long objId = data.ids[index];
			if (data.is2dClassification)
			{
				int ref = data.md.getValueInt(MDLabel.MDL_REF, objId);
				long count = data.md.getValueLong(MDLabel.MDL_CLASS_COUNT, objId);
				return String.format("class %d (%d images)", ref, count);
			}
			else
				return data.md.getValueString(displayLabel.getLabel(), objId);
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * Function to create an image item
	 */
	protected ImageItem createImageItem(int index, int renderLabel, int displayLabel, String key) throws Exception
	{
		String imageFn = getImageFilename(index, renderLabel);
		long objId = data.ids[index];
		ImageItem item = new ImageItem(index);
		boolean useGeo = data.useGeo && renderLabel == MDLabel.MDL_IMAGE;
		ImagePlus imp = getImage(objId, imageFn, thumb_width, thumb_height, useGeo, data.wrap);
		item.setImagePlus(imp);
		return item;
	}

	@Override
	public String getItemKey(int index) throws Exception
	{
		return getItemKey(index, renderLabel.getLabel());

	}

	@Override
	public String getTitle()
	{
		String title = "Metadata: " + (data.getFileName() != null ? Filename.getBaseName(data.getFileName()) : "");
		if (n > 1)
			title += String.format(" %d items", n);
		if (data.hasRenderLabel())
			title += String.format(" (%d x %d)", image_width, image_height);

		return title;
	}// function getTitle

	public String getImageFilename(int index, int label)
	{
		try
		{
			String file = data.getValueFromLabel(index, label);
			String mddir = data.md.getBaseDir();
			file = Filename.findImagePath(file, mddir, true);
			return file;
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		return null;
	}// function getImageFilename

	// @Override
	// public String getImageFilenameAt(int row, int col) {
	// return data.isImageFile(renderLabel) ? data.getValueFromCol(row,
	// renderLabel) : null;
	// }

	@Override
	public boolean handleDoubleClick(int row, int col)
	{
		try
		{
			if (data.isImageFile(renderLabel))
			{
                            int index = getIndex(row, col);
                            String file = getImageFilename(index, renderLabel.getLabel());
                            openXmippImageWindow(file);
                            return true;
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		return false;
	}// function handleDoubleClick
        
        protected void openXmippImageWindow(String file)
        {
                boolean allowsGeometry = data.md.containsGeometryInfo();
                ImagePlusLoader loader = new ImagePlusLoader(file, allowsGeometry, data.useGeo, data.wrap);
                
                if (getNormalized())
                    loader.setNormalize(normalize_min, normalize_max);
                ImagesWindowFactory.openXmippImageWindow(data.window, loader, loader.allowsPoll());
        }

	@Override
	protected double[] getMinAndMax()
	{
		try
		{
			return data.md.getStatistics(false);
		}
		catch (Exception ex)
		{
			DEBUG.printException(ex);
		}
		return null;
	}// function getMinAndMax

	/** Change the use of geometry info */
	public void setUseGeometry(boolean geo, boolean wrap)
	{
		if (!geo)
			wrap = false;
		boolean changed = data.useGeo != geo || data.wrap != wrap;
		data.useGeo = geo;
		data.wrap = wrap;
		if (changed)
			fireTableDataChanged();
	}// function setUseGeometry

	@Override
	public ImagePlusLoader getImageLoader()
	{
		try
		{
                    ImagePlus imp = XmippImageConverter.readMetadataToImagePlus(renderLabel.getLabel(), data.md, data.useGeo, data.wrap);
                    return new ImagePlusLoader(imp);
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		return null;
	}// function getImageLoader

	@Override
	public boolean handleRightClick(int row, int col, XmippPopupMenuCreator xpopup)
	{
		return true;
	}// function handleRightClick

	// Extension of the ImagePlusLoader, read an image from a Metadata row
	


}
