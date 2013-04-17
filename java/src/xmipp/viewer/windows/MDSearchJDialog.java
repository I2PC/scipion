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

import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTable;
import javax.swing.JTextField;

import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippResource;
import xmipp.utils.XmippWindowUtil;

import xmipp.viewer.models.MetadataTableModel;



/* This class will serve for find and replace values
 * in a MetadataGallery
 */
public class MDSearchJDialog extends XmippDialog {
	private static final long serialVersionUID = 1L;
	protected GalleryJFrame parent;
	protected MetaData md;
	protected JTable table;
	protected long[] ids;
	protected int selected_index = -1;
	protected int label;
	protected int[] labels;

	protected GridBagConstraints gbc = new GridBagConstraints();
	protected JPanel mainPanel;
	protected JComboBox jcbLabel, jcbFillMode;

	private JTextField jtFind;
	private JTextField jtReplace;
	private JButton btnReplace;
	private JButton btnFind;
	private JButton btnReplaceFind;
	private JButton btnReplaceAll;
	private JRadioButton jrbForward;
	private JCheckBox jchCase;

	/**
	 * Constructor to add a new label the labels present in the metadata are
	 * passed
	 * */
	public MDSearchJDialog(JFrame parent, JTable table, MetaData md) {
		super(parent, "Find and Replace", true);

		this.parent = (GalleryJFrame) parent;
		this.md = md;
		this.table = table;
		this.ids = md.findObjects();
		this.labels = md.getActiveLabels();
		this.label = labels[0];
		btnOkText = "Close";
		btnCancelDisplay = false;
		initComponents();
	}// constructor MDSearchJDialog

	/**
	 * Utility function to create a pair of label and other component If the
	 * component is null, a textbox will be created and added to the dict
	 */
	protected void addPair(JPanel panel, String text, Component c, int col,
			int row) {
		JLabel label = new JLabel(text);
		gbc.anchor = GridBagConstraints.EAST;
		panel.add(label, XmippWindowUtil.getConstraints(gbc, col, row));
		gbc.anchor = GridBagConstraints.WEST;
		panel.add(c, XmippWindowUtil.getConstraints(gbc, col + 1, row));
	}

	@Override
	protected void createContent(JPanel panel) {
		setResizable(false);
		mainPanel = new JPanel(new GridBagLayout());
		gbc.insets = new Insets(5, 5, 0, 0);

		jtFind = new JTextField(15);
		jtReplace = new JTextField(15);

		int row = 0;
		addPair(mainPanel, "Find", jtFind, 0, row++);
		addPair(mainPanel, "Replace", jtReplace, 0, row++);

		jtFind.addActionListener(this);

		// setPreferredSize(new Dimension(350, 180));
		setIconImage(XmippResource.getIcon("search.gif").getImage());

		JPanel optPanel = new JPanel(new GridBagLayout());
		mainPanel.add(optPanel,
				XmippWindowUtil.getConstraints(gbc, 0, row++, 2));
		optPanel.setBorder(BorderFactory.createTitledBorder("Options"));

		// Create labels combo
		jcbLabel = new JComboBox();
		for (int l : labels) {
			jcbLabel.addItem(MetaData.getLabelName(l));
		}
		jcbLabel.addActionListener(this);
		addPair(optPanel, "Label", jcbLabel, 0, 0);
		jchCase = new JCheckBox();
		addPair(optPanel, "Case sensitive", jchCase, 0, 1);
		ButtonGroup btnGroup = new ButtonGroup();
		jrbForward = new JRadioButton("Forward");
		jrbForward.setSelected(true);
		btnGroup.add(jrbForward);
		JRadioButton jrb = new JRadioButton("Backward");
		btnGroup.add(jrb);
		addPair(optPanel, "Direction", jrbForward, 0, 2);
		addPair(optPanel, "", jrb, 0, 3);

		// Add buttons
		JPanel btnPanel = new JPanel(new GridBagLayout());

		btnFind = XmippWindowUtil.getTextButton("Find", this);
		btnReplace = XmippWindowUtil.getTextButton("Replace", this);
		btnReplaceFind = XmippWindowUtil.getTextButton("Replace/Find", this);
		btnReplaceAll = XmippWindowUtil.getTextButton("Replace All", this);

		btnPanel.add(btnFind, XmippWindowUtil.getConstraints(gbc, 0, 0, 2));
		btnPanel.add(btnReplaceFind,
				XmippWindowUtil.getConstraints(gbc, 2, 0, 2));
		btnPanel.add(btnReplace, XmippWindowUtil.getConstraints(gbc, 0, 1, 2));
		btnPanel.add(btnReplaceAll,
				XmippWindowUtil.getConstraints(gbc, 2, 1, 2));
		gbc.anchor = GridBagConstraints.EAST;
		mainPanel.add(btnPanel,
				XmippWindowUtil.getConstraints(gbc, 0, row++, 2));

		panel.add(mainPanel);
	}// function initComponents

	@Override
	public void handleActionPerformed(ActionEvent evt) {
		Object o = evt.getSource();
		if (o == jtFind || o == btnFind)
			find();
		else if (o == jcbLabel)
			label = labels[jcbLabel.getSelectedIndex()];
		else if (o == btnReplace)
			replace(selected_index);
		else if (o == btnReplaceFind) {
			replace(selected_index);
			find();
		} else if (o == btnReplaceAll)
			replaceAll();

	}// function actionPerformed

	private int nextIndex(int index, int direction) {
		int next = 0;
		if (index >= 0) {
			next = (index + direction);
			if (next < 0)
				next = ids.length - 1;
			else if (next >= ids.length)
				next = 0;
		}
		return next;
	}

	private void find() {
		try {
			int direction = jrbForward.isSelected() ? 1 : -1;
			String findStr = jtFind.getText();
			if (findStr.length() == 0)
				return;

			int index = nextIndex(selected_index, direction);
			selected_index = -1;

			for (int i = 0; i < ids.length && selected_index < 0; ++i) {
				if (match(findStr, index))
					selected_index = index;
				else
					index = nextIndex(index, direction);
			}
			if (selected_index > -1)
				parent.selectItem(selected_index, 0);
		} catch (Exception e) {
			showException(e);
		}
	}

	private boolean match(String findStr, int index) {
		boolean caseSensitive = jchCase.isSelected();
		String value = md.getValueString(label, ids[index]).trim();
		String str = findStr.trim();
		if (!caseSensitive) {
			value = value.toLowerCase();
			str = findStr.toLowerCase();
		}
		return (value.indexOf(str) >= 0);
	}

	private void replace(int index) {
		String replaceStr = jtReplace.getText().trim();
		if (index > -1 && replaceStr.length() > 0) {
			String value = md.getValueString(label, ids[index]);
			String newValue = value.replace(jtFind.getText().trim(), replaceStr);
			// table.setValueAt(newValue, index, jcbLabel.getSelectedIndex());
			md.setValueString(label, newValue, ids[index]);
			((MetadataTableModel) table.getModel()).fireTableRowsUpdated(index,
					index);
		}
	}

	private void replaceAll() {
		String findStr = jtFind.getText().trim();
		for (int i = 0; i < ids.length; ++i)
			if (match(findStr, i))
				replace(i);
	}

}// class AddFillLabelsJDialog
