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

package xmipp.ij.commons;

import java.io.File;

import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import ij.ImagePlus;

public class ImagePlusLoader
{

	protected String fileName = null;
	protected boolean allowsPoll;
	protected boolean allowsGeometry = false;
	protected boolean useGeometry;
	protected ImagePlus imp;
	protected ImageGeneric ig;
	protected long modified;
	protected boolean wrap;
	
	public ImagePlusLoader()
	{
		this.imp = null;
	}

	public ImagePlusLoader(ImagePlus imp)
	{
		this.imp = imp;
		allowsPoll = existsFile();
	}

	public ImagePlusLoader(String fileName)
	{
		this.fileName = fileName;
		this.modified = new File(fileName).lastModified();
		allowsPoll = existsFile();
	}

	public ImagePlusLoader(ImageGeneric ig)
	{
		this.ig = ig;
		allowsPoll = existsFile();
		if(existsFile())
			fileName = ig.getFilename();
	}

	public ImagePlus getImagePlus()
	{
		if (imp == null)
			imp = loadImagePlus();
		return imp;
	}

	public ImagePlus loadImagePlus()
	{
		imp = null;
		try
		{
			if (fileName != null && Filename.exists(fileName) && (hasChanged() || imp == null))
				imp = loadImage();
			else if (ig != null)
				imp = XmippImageConverter.readToImagePlus(ig);
			return imp;
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		return imp;
	}

	protected ImagePlus loadImage() throws Exception
	{
		ig = new ImageGeneric(fileName);
		long select_image = Filename.getNimage(fileName);
		imp = XmippImageConverter.readToImagePlus(ig, ig.getXDim(), ig.getYDim(), select_image);
		return imp;
		
	}

	public boolean hasChanged()
	{
		return new File(fileName).lastModified() > modified;
	}

	public String getFileName()
	{
		return fileName;
	}

	public boolean allowsPoll()
	{
		return allowsPoll;
	}

	public boolean allowsGeometry()
	{
		return allowsGeometry;
	}

	public void setAllowsGeometry(boolean useGeometry)
	{
		this.useGeometry = useGeometry;
	}

	public boolean getUseGeometry()
	{
		return useGeometry;
	}

	public void setUseGeometry(boolean value)
	{
		useGeometry = value;
		
	}

	public void setWrap(boolean value)
	{
		wrap = value;

	}

	public boolean isWrap()
	{
		return wrap;
	}


	public boolean isVolume()
	{
		try
		{
			if (ig == null)
				return false;
			return ig.isVolume();

		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}
	}
	
	public boolean existsFile()
	{
		String file = null;
		if(fileName != null && !fileName.equals(""))
			file = fileName;
		else if(imp != null && imp.getOriginalFileInfo() != null)
			file = imp.getOriginalFileInfo().directory + File.separator + imp.getOriginalFileInfo().fileName;
		else if (ig!= null && ig.getFilename()!= null)
			file = ig.getFilename();
		if(file == null)
			return false;
		if(!new File(file).exists())
			return false;
		return true;
	}

}
