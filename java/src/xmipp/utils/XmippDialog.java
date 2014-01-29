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
import java.awt.Dialog;
import java.awt.FlowLayout;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

/**
 * This class serve as a simple JDialog extension. The Ok and Cancel buttons are
 * added and the dispose of the dialog is set. After constructing a dialog the
 * show function should be called. It will return true if OK is pressed and
 * false if Cancel.
 */
public class XmippDialog extends JDialog implements ActionListener {
	private static final long serialVersionUID = 1L;
	protected Window parent;
	protected boolean result;
	protected JButton btnCancel;
	protected JButton btnOk;
	protected String btnOkText = "Ok";
	protected String btnCancelText = "Cancel";
	protected boolean btnCancelDisplay = true;
	protected JPanel panelBtn;
	protected boolean disposeOnClose = true;
	protected boolean closeOnAction = true;
        
        

	public XmippDialog(JFrame parent, String title, boolean modal) {
		super(parent, title, modal);
		this.parent = parent;
	}

	public XmippDialog(Dialog parent, String title, boolean modal) {
		super(parent, title, modal);
		this.parent = parent;
	}

	/**
	 * this is the general method to init the components. It should be called
	 * from every subclass constructor after some settings of values
	 */
	protected void initComponents() {
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		Container content = getContentPane();
		JPanel panel = new JPanel();
		createContent(panel);
		content.add(panel, BorderLayout.CENTER);
		createButtons();
		content.add(panelBtn, BorderLayout.PAGE_END);
		pack();
		if (parent != null)
			XmippWindowUtil.centerWindows(this, parent);
	}

	protected JButton addButton(String text) {
		JButton btn = XmippWindowUtil.getTextButton(text, this);
		panelBtn.add(btn);
		return btn;
	}

	protected void createButtons() {
		// Create panel for Ok and Cancel buttons
		panelBtn = new JPanel(new FlowLayout(FlowLayout.TRAILING));
		btnOk = addButton(btnOkText);
		if (btnCancelDisplay)
			btnCancel = addButton(btnCancelText);
	}

	/** Function to display the Dialog and return the result state */
	public boolean showDialog() {
		setVisible(true);
		return result;
	}

	/**
	 * This function should be redefined to add the content to the dialog, the
	 * container panel is passed
	 */
	protected void createContent(JPanel panel) {

	}

	/**
	 * This function should be overrided to handle action events.
	 */
	public void handleActionPerformed(ActionEvent evt) {

	}

	/**
	 * Function called before close after ok pressed
	 */
	public void handleOk() {
	}

	/**
	 * Function called before close after cancel pressed
	 */
	public void handleCancel() {
	}
        
        

	@Override
	public void actionPerformed(ActionEvent evt) {
		Object obj = evt.getSource();
		if (obj == btnOk) {
			close(true);
		} else if (obj == btnCancel) {
			close(false);
		} else
			handleActionPerformed(evt);
	}

	/** Close function to hide the Dialog and dispose resources */
	protected void close(boolean result) {
		if (result)
			handleOk();
		else
			handleCancel();
		
		if (closeOnAction) {
			this.result = result;
			setVisible(false);
			if (disposeOnClose)
				dispose();
		}
	}// function close

	/** Change the default dispose on close behavior */
	public void setDisposeOnClose(boolean value) {
		disposeOnClose = value;
	}

	/** Set to false if don't want to close on close */
	public void setCloseOnOk(boolean value) {
		closeOnAction = value;
	}

	/** Some methods to show dialog using the current as parent */
	public boolean showInfo(String message) {
		XmippMessageDialog dlg = new XmippMessageDialog(this, "INFO", message,
				"info.gif");
		return dlg.showDialog();
	}

	public boolean showWarning(String message) {
		XmippMessageDialog dlg = new XmippMessageDialog(this, "WARNING",
				message, "warning.gif");
		return dlg.showDialog();
	}

	public boolean showError(String message) {
		XmippMessageDialog dlg = new XmippMessageDialog(this, "ERROR", message,
				"error.gif");
		return dlg.showDialog();
	}

	public boolean showException(Exception e) {
		XmippMessageDialog dlg = new XmippMessageDialog(this, "ERROR",
				e.getMessage(), "error.gif");
		// TODO: Integrate the stack trace into the dialog
		e.printStackTrace();
		return dlg.showDialog();
	}

	/** Some static methods to display some message dialogs */
	public static boolean showInfo(JFrame parent, String message) {
		XmippMessageDialog dlg = new XmippMessageDialog(parent, "INFO",
				message, "info.gif");
		dlg.btnCancel.setVisible(false);
		return dlg.showDialog();
	}

	public static boolean showWarning(JFrame parent, String message) {
		XmippMessageDialog dlg = new XmippMessageDialog(parent, "WARNING",
				message, "warning.gif");
		return dlg.showDialog();
	}

	public static boolean showError(JFrame parent, String message) {
		XmippMessageDialog dlg = new XmippMessageDialog(parent, "ERROR",
				message, "error.gif");
		return dlg.showDialog();
	}

	public static boolean showException(JFrame parent, Exception e) {
		XmippMessageDialog dlg = new XmippMessageDialog(parent, "ERROR",
				e.getMessage(), "error.gif");
		// TODO: Integrate the stack trace into the dialog
		e.printStackTrace();
		return dlg.showDialog();
	}

	public static boolean showQuestion(JFrame parent, String message) {
		XmippMessageDialog dlg = new XmippQuestionDialog(parent, message);
		return dlg.showDialog();
	}

	public static boolean showQuestion(JDialog parent, String message) {
		XmippMessageDialog dlg = new XmippQuestionDialog(parent, message);
		return dlg.showDialog();
	}
        
        public static int showQuestionYesNoCancel(JFrame parent, String message) {
		XmippQuestionDialog dlg = new XmippQuestionDialog(parent, message);
		dlg.setVisible(true);
                return dlg.getOption();
	}

	public static String removeColors(String message) {
		if (message == null)
			return null;
		String redPrefix = String.format("%c[%d;%dm", 0x1B, 1, 31);
		String redSuffix = String.format("%c[0m", 0x1B);
		return message.replace(redPrefix, "").replace(redSuffix, "");
	}
        
        

}// class XmippDialog

