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
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;

import javax.swing.BorderFactory;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.models.MetadataGalleryTableModel;


public class AddFillLabelsJDialog extends XmippDialog {
	private static final long serialVersionUID = 1L;
	protected MetadataGalleryTableModel gallery;

	protected boolean fillMode = false;
	protected GridBagConstraints gbc = new GridBagConstraints();
	protected int label;
	protected ArrayList<ColumnInfo> labels = null;
	protected Hashtable<String, JLabel> dictLabels = new Hashtable<String, JLabel>();
	protected Hashtable<String, JTextField> dictTexts = new Hashtable<String, JTextField>();
	protected Hashtable<String, JPanel> dictPanels = new Hashtable<String, JPanel>();
	protected JPanel mainPanel;
	protected JComboBox jcbLabel, jcbFillMode;

	/**
	 * Constructor to add a new label the labels present in the metadata are
	 * passed
	 * */
	public AddFillLabelsJDialog(GalleryJFrame parent,
			ArrayList<ColumnInfo> labels) {
		super(parent, "Add new label", true);
		this.gallery = (MetadataGalleryTableModel) parent.gallery;
		this.labels = labels;
		initComponents();
	}// constructor AddFillLabelsJDialog

	/** Constructor to fill an existing label */
	public AddFillLabelsJDialog(GalleryJFrame parent, int label) {
		super(parent, "Fill label values", true);
		this.label = label;
		initComponents();
	}// constructor AddFillLabelsJDialog

	/** Return if in fill mode */
	public boolean isFillMode() {
		return this.labels == null;
	}

	/**
	 * Utility function to create a pair of label and other component If the
	 * component is null, a textbox will be created and added to the dict
	 */
	protected void addPair(JPanel panel, String text, int row, Component c) {
		JLabel label = new JLabel(text);
		dictLabels.put(text, label);
		gbc.anchor = GridBagConstraints.EAST;
		panel.add(label, XmippWindowUtil.getConstraints(gbc, 0, row));
		if (c == null) {
			c = new JTextField(10);
			dictTexts.put(text, (JTextField) c);
		}
		gbc.anchor = GridBagConstraints.WEST;
		panel.add(c, XmippWindowUtil.getConstraints(gbc, 1, row));
	}

	/** Create a combobox with possible new labels */
	private JComboBox createLabelsCombo() {
		jcbLabel = new JComboBox();
		boolean found = false;
		for (int label = MDLabel.MDL_OBJID + 1; label < MDLabel.MDL_LAST_LABEL; ++label) {
			found = false;
			for (ColumnInfo ci : this.labels)
				if (ci.getLabel() == label) {
					found = true;
					break;
				}
			if (!found) {
				String name = MetaData.getLabelName(label);
				if (!name.startsWith("emx_"))
					jcbLabel.addItem(name);
			}
		}
		jcbLabel.addActionListener(this);
		return jcbLabel;
	}

	/** Create panel with values pairs needed */
	private JPanel addPairsPanel(String option, int panelRow, String... values) {
		JPanel panel = new JPanel(new GridBagLayout());
		int row = 0;
		for (String value : values) {
			addPair(panel, value, row++, null);
		}
		mainPanel.add(panel,
				XmippWindowUtil.getConstraints(gbc, 0, panelRow, 2));
		dictPanels.put(option, panel);
		panel.setBorder(BorderFactory.createTitledBorder("Params"));
		panel.setVisible(false);
		return panel;
	}

	@Override
	protected void createContent(JPanel panel) {
		setResizable(false);
		mainPanel = new JPanel(new GridBagLayout());
		gbc.insets = new Insets(5, 5, 0, 0);

		Component c = null;
		if (isFillMode()) {
			JTextField jtext = new JTextField(MetaData.getLabelName(label));
			jtext.setEditable(false);
			c = jtext;
			setPreferredSize(new Dimension(250, 200));
		} else {
			c = createLabelsCombo();
			setPreferredSize(new Dimension(350, 200));
		}
		addPair(mainPanel, "Label", 0, c);
		String[] opts = { MetaData.FILL_CONSTANT, MetaData.FILL_LINEAR, MetaData.FILL_RAND_UNIFORM,
				MetaData.FILL_RAND_GAUSSIAN };
		jcbFillMode = new JComboBox(opts);
		jcbFillMode.addActionListener(this);
		addPair(mainPanel, "Fill mode", 1, jcbFillMode);
		addPairsPanel(opts[0], 2, "value").setVisible(true);
		addPairsPanel(opts[1], 3, "start", "step");
		addPairsPanel(opts[2], 4, "min", "max");
		addPairsPanel(opts[3], 5, "mean", "std");

		panel.add(mainPanel);

	}// function initComponents

	/** some helper functions */
	public String getFillMode() {
		return (String) jcbFillMode.getSelectedItem();
	}

	/** Return the current working label */
	public int getLabel() {
		if (isFillMode())
			return label;
		else {
			try {
				return MetaData.str2Label((String) jcbLabel.getSelectedItem());
			} catch (Exception e) {
				showException(e);
			}
			return MDLabel.MDL_UNDEFINED;
		}
	}

	public String getValue(String key) {
		JTextField jtf = dictTexts.get(key);
		return jtf.getText();
	}
	
	public String[] getValues(){
		String[] values = null;
		String mode = getFillMode();
		if (mode.equalsIgnoreCase(MetaData.FILL_CONSTANT))
			values = new String[] { getValue("value") };
		else if (mode.equalsIgnoreCase(MetaData.FILL_LINEAR))
			values = new String[] { getValue("start"), getValue("step") };
		else if (mode.equalsIgnoreCase(MetaData.FILL_RAND_UNIFORM))
			values = new String[] { getValue("min"), getValue("max") };
		else if (mode.equalsIgnoreCase(MetaData.FILL_RAND_GAUSSIAN))
			values = new String[] { getValue("mean"), getValue("std") };
		
		return values;
	}

	@Override
	public void handleActionPerformed(ActionEvent evt) {
		JComboBox jcb = (JComboBox) evt.getSource();
		if (jcb == jcbFillMode) {
			String fillMode = getFillMode();
			for (Enumeration<String> keys = dictPanels.keys(); keys
					.hasMoreElements();) {
				String keyMode = keys.nextElement();
				dictPanels.get(keyMode).setVisible(fillMode.equals(keyMode));
			}
		}
	}// function actionPerformed

	/** Perform the selected action */
//	@Override
//	public void handleOk() {
//		try {
//			gallery.fillLabel(getLabel(), getFillMode(), dictTexts.get("value").getText());
//		} catch (Exception e) {
//			XmippDialog.showException(parent, e);
//		}
//	}

}// class AddFillLabelsJDialog
