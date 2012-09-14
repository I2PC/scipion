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
import java.awt.Color;
import java.awt.Dialog;
import java.awt.FlowLayout;

import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

/** Special case of XmippDialog to display message (info, error, warning) */
public class XmippMessageDialog extends XmippDialog {
	private static final long serialVersionUID = 1L;
	String message;
	String iconPath;

	public XmippMessageDialog(JFrame parent, String title){
		super(parent, title, true);
	}
	
	public XmippMessageDialog(Dialog parent, String title){
		super(parent, title, true);
	}
	
	public XmippMessageDialog(JFrame parent, String title, String message,
			String iconPath) {
		super(parent, title, true);
		init(removeColors(message), iconPath);
	}
	
	public XmippMessageDialog(Dialog parent, String title, String message,
			String iconPath) {
		super(parent, title, true);
		init(message, iconPath);
	}
	
	/** Perform some initializations */
	protected void init(String message, String iconPath){
		this.message = message;
		this.iconPath = iconPath;
		initComponents();
	}

	@Override
	protected void createContent(JPanel panel) {
		panel.setLayout(new BorderLayout());
		JPanel iconPanel = new JPanel();
		iconPanel.add(XmippWindowUtil.getIconLabel(iconPath));
		iconPanel.setBackground(Color.white);
		panel.add(iconPanel, BorderLayout.LINE_START);
		JTextArea text = new JTextArea(message, 5, 40);
		text.setEditable(false);
		text.setBackground(Color.white);
		text.setBorder(BorderFactory.createEmptyBorder());
		text.setLineWrap(true);
		text.setWrapStyleWord(true);
		// text
		JScrollPane scrollPane = new JScrollPane(text);
		scrollPane.setBorder(BorderFactory.createEmptyBorder());
		scrollPane.setBackground(Color.white);
		panel.add(scrollPane, BorderLayout.CENTER);
	}

	@Override
	protected void createButtons() {
		super.createButtons();
		panelBtn.setLayout(new FlowLayout());
	}
}// class XmippMessageDialog