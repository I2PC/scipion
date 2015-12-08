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

import java.util.ArrayList;
import java.util.List;

import xmipp.ij.commons.ImagePlusLoader;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.utils.DEBUG;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.viewer.ImageDimension;
import xmipp.viewer.windows.ImagesWindowFactory;

//Gallery mode!!!!!!!!!!
public class MetadataGalleryTableModel extends ImageGalleryTableModel
{

	private static final long serialVersionUID = 1L;
	// Label to be rendered
	protected ColumnInfo renderLabel;
	// Also store the visible ones to fast access
	public ArrayList<ColumnInfo> visibleLabels;

	public MetadataGalleryTableModel(GalleryData data, boolean[] selection) throws Exception
	{
		super(data, selection);
		data.normalize = false;
	}

	/** Update the columns display information */
	public void updateColumnInfo(List<ColumnInfo> newInfo)
	{
		int n = newInfo.size();
		boolean changed = false;
		data.renderImages = false;
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
						data.renderImages = true;
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
	
	public ImagePlus getImage(long id, String imagepath, ColumnInfo renderLabel)
	{
		ImagePlus imp = null;
		if (imagepath != null && Filename.exists(imagepath))
		{
			try
			{
                ImagePlusLoader loader;
                if(Filename.isStackOrVolume(imagepath))
                    loader = new ImagePlusLoader(imagepath, null, null, data.useGeo, data.wrap, data.inverty, ImageGeneric.MID_SLICE);
                else  
                    loader = new ImagePlusLoader(imagepath, data.useGeo, data.wrap, data.inverty);
                loader.setGeometry(data.getGeometry(id, renderLabel));
                loader.setInvertY(data.inverty);
                loader.setDimension(thumb_width, thumb_height);
                return loader.getImagePlus();
			}
			catch (Exception ex)
			{
				imp = null;
			}
		}
		return imp;
	}
        
        
    protected void openXmippImageWindow(int index, ColumnInfo ci)
    {
        String file = getImageFilename(index, ci.label);
        if(file == null)
            throw new IllegalArgumentException(XmippMessage.getPathNotExistsMsg(data.getValueFromLabel(index, ci.label)));
        ImagePlusLoader loader = new ImagePlusLoader(file, data.useGeo, data.wrap, data.inverty);
        if(data.containsGeometryInfo())
            loader.setGeometry(data.getGeometry(data.ids[index], ci));
        if (data.getNormalized())
            loader.setNormalize(normalize_min, normalize_max);
        ImagesWindowFactory.openXmippImageWindow(loader, data.parameters);
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
        int width = 50, height = 50;
		if (data.hasRenderLabel() && data.renderImages) 
		{
			renderLabel = data.ciFirstRender;
			// if (renderLabels) {
			String imageFn = data.getSampleImage(renderLabel);
            if(imageFn != null)
            {
                	ImagePlusLoader loader;
                	if(Filename.isStackOrVolume(imageFn))
                        loader = new ImagePlusLoader(imageFn, null, null, false, false, false, ImageGeneric.MID_SLICE);
                    else  
                        loader = new ImagePlusLoader(imageFn);
                    ImagePlus image = loader.getImagePlus();
                    if(image != null && image.getWidth() > 0)
                    {
	                    width = image.getWidth(); 
	                    height = image.getHeight();
	                    dim = new ImageDimension(width, height);
                    }
            }
            
		}
		
		if (dim == null)
			dim = new ImageDimension(width);
		// Zdim will always be used as number of elements to display
		dim.setZDim(data.ids.length);
		n = dim.getZDim();
		image_width = dim.getXDim();
		image_height = dim.getYDim();
		return dim;
	}

	@Override
	protected ImageItem createItem(int index, String key, ColumnInfo ci) throws Exception
	{
		return createImageItem(index, ci);
	}
	
	

	public String getLabel(int row, int col)
	{
		try
		{
			int index = getIndex(row, col);
			long objId = data.ids[index];
            return data.getDisplayLabel(objId);
                        
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
	protected ImageItem createImageItem(int index, ColumnInfo renderLabel) throws Exception
	{
		String imageFn = getImageFilename(index, renderLabel.label);
		long objId = data.ids[index];
		ImagePlus imp = getImage(objId, imageFn, renderLabel);
		ImageItem item = new ImageItem(index, imp);
		return item;
	}

	@Override
	public String getItemKey(int index) throws Exception
	{
		return getItemKey(index, renderLabel.label);

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
                openXmippImageWindow(index, getColumn(row, col));
                return true;
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		return false;
	}// function handleDoubleClick
        
        
	@Override
	protected double[] getMinAndMax()
	{
		try
		{
			return data.md.getStatistics(false, data.ciFirstRender.label);
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
            ImagePlus imp = XmippImageConverter.readMetadataToImagePlus(renderLabel.label, data.md, data.useGeo, data.wrap, data.inverty);
            return new ImagePlusLoader(imp, data.inverty);
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

    @Override
    public boolean showLabels() {
        return true;
    }
	
    
	
        

}
