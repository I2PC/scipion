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

package xmipp.utils;

import java.io.File;
import javax.swing.*;
import javax.swing.filechooser.*;

import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;

/* This class implements the abstract FileView to personalize
 * the default icons and file info displayed on XmippFileChooser */
public class XmippFileView extends FileView {

	public String getName(File f) {
		return null; // let the L&F FileView figure this out
	}

	public String getDescription(File f) {
		try {
			if (!f.isDirectory()) {
				String filename = f.getPath();
				if (Filename.isSingleImageExt(filename)
						|| Filename.isStackExt(filename)
						|| Filename.isVolumeExt(filename)) {
					ImageGeneric image = new ImageGeneric(filename);
					String desc = String.format("<Image dimensions: %d x %d",
							image.getXDim(), image.getYDim());
					return desc;
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null; // let the L&F FileView figure this out
	}

	public Boolean isTraversable(File f) {
		return null; // let the L&F FileView figure this out
	}

	public String getTypeDescription(File f) {
		String type = null;
		return type;
	}

	public Icon getIcon(File f) {
		Icon icon = null;
		try {
			String filename = f.getPath();
			// if (Filename.exists(filename))
			{
				String iconString = "generic_file.gif";
				if (f.isDirectory())
					iconString = "folderopen.gif";
				else if (Filename.isSingleImageExt(filename))
					iconString = "image.gif";
				else if (Filename.isMetadataExt(filename)) 
					iconString = "md.gif";
				else if (Filename.isVolumeExt(filename))
					iconString = "vol.gif";
				else if (Filename.isStackExt(filename))
					iconString = "stack.gif";
				else if (filename.endsWith(Filename.EXT_ERR))
					iconString = "err.gif";
				else if (filename.endsWith(Filename.EXT_OUT))
					iconString = "out.gif";
				else if (filename.endsWith(Filename.EXT_LOG))
					iconString = "log.gif";
				else if (filename.endsWith(Filename.EXT_PY))
					iconString = "python_file.gif";
				return XmippResource.getIcon(iconString);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return icon;
	}
}
