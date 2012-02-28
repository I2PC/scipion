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

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JPanel;

/** This class serve as a simple JDialog extension.
 * The Ok and Cancel buttons are added and the dispose
 * of the dialog is set. After constructing a dialog
 * the show function should be called. It will return
 * true if OK is pressed and false if Cancel.
 */
public class XmippDialog extends JDialog  implements ActionListener {
	private static final long serialVersionUID = 1L;
	protected JFrame parent;
	protected boolean result;
	protected JButton btnCancel;
	protected JButton btnOk;
	protected String okText = "Ok";
	protected String cancelText = "Cancel";
	
	public XmippDialog(JFrame parent, String title, boolean modal){
		super(parent, title, modal);
		this.parent = parent;
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle(title);
		Container content = getContentPane();
		JPanel panel = new JPanel();
		createContent(panel);
		content.add(panel, BorderLayout.CENTER);
		//Create panel for Ok and Cancel buttons
		JPanel panelBtn = new JPanel(new FlowLayout(FlowLayout.TRAILING));
		btnOk = WindowUtil.getTextButton(okText, this);
		panelBtn.add(btnOk);
		btnCancel = WindowUtil.getTextButton(cancelText, this);
		panelBtn.add(btnCancel);
		content.add(panelBtn, BorderLayout.PAGE_END);
		pack();
		WindowUtil.centerWindows(this, parent);
	}
	
	/** Function to display the Dialog and return the result state */
	public boolean showDialog() {
		setVisible(true);
		return result;
	}
	
	/** This function should be redefined to add the content
	 * to the dialog, the container panel is passed */
	protected void createContent(JPanel panel){
		
	}
	
	/** This function should be overrided to handle action
	 * events.
	 */
	public void handleActionPerformed(ActionEvent evt){
		
	}

	@Override
	public void actionPerformed(ActionEvent evt) {
		Object obj = evt.getSource();
		if (obj == btnOk)
			close(true);
		else if (obj == btnCancel) {
			close(false);
		}
		else 
			handleActionPerformed(evt);		
	}
	
	/** Close function to hide the Dialog and dispose resources */
	protected void close(boolean result) {
		this.result = result;
		setVisible(false);
		dispose();
	}// function close
	
}
