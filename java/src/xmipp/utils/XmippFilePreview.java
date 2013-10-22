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

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;

import javax.swing.ImageIcon;
import javax.swing.JComponent;

import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MetaData;

/* ImagePreview.java by FileChooserDemo2.java. */
public class XmippFilePreview extends JComponent implements
		PropertyChangeListener {
	ImageIcon thumbnail = null;
	File file = null;

	public XmippFilePreview(XmippFileChooser fc) {
		setPreferredSize(new Dimension(128, 128));
		fc.addPropertyChangeListener(this);
	}

	private void loadImage(String filename) throws Exception{
		try{
		ImageGeneric image;
		image = new ImageGeneric(filename);
		image.read(128, 128, ImageGeneric.MID_SLICE, ImageGeneric.FIRST_IMAGE);
		image.convert2Datatype(ImageGeneric.UChar);
		BufferedImage bimg = new BufferedImage(image.getXDim(), image.getYDim(),
				BufferedImage.TYPE_BYTE_GRAY);
		bimg.getRaster().setDataElements(0, 0, image.getXDim(), image.getYDim(), 
				image.getArrayByte(ImageGeneric.FIRST_IMAGE, ImageGeneric.FIRST_SLICE));
		thumbnail = new ImageIcon(bimg);
		} catch (Exception e){
			thumbnail = null;
		}
	}
	
	public void loadImage() {
		if (file == null) {
			thumbnail = null;
			return;
		}

		try {
			String filename = file.getPath();
			if (!file.exists() || file.isDirectory()){
				thumbnail = null;
				return;
			}
			if (Filename.isSingleImageExt(filename) ||
					Filename.isStackExt(filename) ||
					Filename.isVolumeExt(filename))
				loadImage(filename);
			else if (Filename.isMetadataExt(filename)){
				MetaData md = new MetaData(filename);
				int[] labels = md.getActiveLabels();
				for (int l: labels)
					if (MetaData.isImage(l))
						loadImage(md.getValueString(l, md.firstObject()));
				md.destroy();
			}			
				
           
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void propertyChange(PropertyChangeEvent e) {
		boolean update = false;
		String prop = e.getPropertyName();

		// If the directory changed, don't show an image.
		if (XmippFileChooser.DIRECTORY_CHANGED_PROPERTY.equals(prop)) {
			file = null;
			update = true;

			// If a file became selected, find out which one.
		} else if (XmippFileChooser.SELECTED_FILE_CHANGED_PROPERTY.equals(prop)) {
			file = (File) e.getNewValue();
			update = true;
		}

		// Update the preview accordingly.
		if (update) {
			thumbnail = null;
			if (isShowing()) {
				loadImage();
				repaint();
			}
		}
	}

	protected void paintComponent(Graphics g) {
		if (thumbnail == null) {
			loadImage();
		}
		if (thumbnail != null) {
			int x = getWidth() / 2 - thumbnail.getIconWidth() / 2;
			int y = getHeight() / 2 - thumbnail.getIconHeight() / 2;

			if (y < 0) {
				y = 0;
			}

			if (x < 5) {
				x = 5;
			}
			thumbnail.paintIcon(this, g, x, y);
		}
	}
}
