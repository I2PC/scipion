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

package xmipp.viewer.windows;

import javax.swing.JPanel;
import javax.swing.JLabel;
import javax.swing.JTextField;
import xmipp.jni.MetaData;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.models.ImageGalleryTableModel;

public class SaveImagesJDialog extends SaveJDialog {
	private static final long serialVersionUID = 1L;
	ImageGalleryTableModel gallery;
	JPanel panelButtons;
	JTextField textField;
	MetaData imagesMd;

	public SaveImagesJDialog(GalleryJFrame parent, String proposedFileName) {
		super(parent, proposedFileName, true);
	}// constructor ColumnsJDialog
	
	@Override
	protected void initComponents(){
		this.gallery = ((GalleryJFrame)parent).gallery;
		this.btnOkText = "Save Images";
		super.initComponents();
	}

	@Override
	protected void createContent(JPanel panel) {
		super.createContent(panel);
                
		imagesMd = gallery.data.getClassesImages(gallery.getSelection());
		String text = String.format("<html>You are about to save <font color='red'>%d</font> images from <font color='red'>%d</font> classes.",
				imagesMd.size(), gallery.getSelectionCount());
		panel.add(new JLabel(text), XmippWindowUtil.getConstraints(gbc, 0, 1, 2));
	}// function initComponents
	
	@Override
	protected void fillMdOptions(){
		// Just to avoid adding some metadata options.
		// this dialog will override the existing file
	}

	@Override
	public void handleOk() {
		imagesMd.write(getMdFilename());
		imagesMd.destroy();
	}// function actionPerformed
	
	@Override
	public void handleCancel(){
		imagesMd.destroy();
	}
}// class ColumnsJDialog
