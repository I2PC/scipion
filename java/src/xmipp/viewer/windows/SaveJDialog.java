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

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.io.File;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.models.GalleryData;

public class SaveJDialog extends XmippDialog
{
	protected static final long serialVersionUID = 1L;
	protected JPanel group;
	protected GridBagConstraints gbc;

	protected XmippFileChooser fc;
	protected JCheckBox chbMd;
//	protected JCheckBox chbImg;
	// protected JButton btnBrowseMd;
	// protected JTextField tbMd;
	// protected JTextField tbImg;
	// protected JButton btnBrowseImg;
	protected BrowseField browseMd;
	protected BrowseField browseImg;
	protected JComboBox cbExtension;
	protected JComboBox cbLabel;
//	protected JRadioButton rbStack;
	protected JRadioButton rbIndependent;
	protected JRadioButton rbMdOverride;
	protected JRadioButton rbMdAppend;
//	protected JPanel panelImg;
	protected JPanel panelMd;
	protected GalleryData data;
	private String block;
	private JCheckBox chbDiscard;
	private boolean isselection;

	public SaveJDialog(GalleryJFrame parent, String file, boolean isselection)
	{
		super(parent, "Save", true);
		this.isselection = isselection;
		initComponents();
		setMdFilename(file);
	}// constructor SaveJDialog

	public SaveJDialog(GalleryJFrame parent)
	{
		this(parent, null, false);
	}// constructor SaveJDialog

	@Override
	protected void createContent(JPanel panel)
	{
		this.data = ((GalleryJFrame) parent).getData();
		//setMinimumSize(new Dimension(500, 300));
		setResizable(false);
		panel.setLayout(new GridBagLayout());
		gbc = new GridBagConstraints();
		gbc.ipadx = 5;
		gbc.ipady = 5;
		gbc.weightx = 1.0;
		gbc.weighty = 1.0;
		// gbc.insets = new Insets(5, 5, 5, 5);
		// gbc.weightx = 1.0;
		// gbc.weighty = 1.0;
		gbc.anchor = GridBagConstraints.NORTHWEST;
		fc = new XmippFileChooser();
		fc.setCurrentDirectory(new File(System.getProperty("user.dir")));

		group = new JPanel(new GridBagLayout());
		group.setBorder(BorderFactory.createTitledBorder("Save options"));

		// Checkbox to save metadata
		// chbMd = new JCheckBox("Save metadata", true);
		// chbMd.addActionListener(this);
		group.add(new JLabel("Metadata path:"), XmippWindowUtil.getConstraints(gbc, 0, 0));

		createMdOptions();
		group.add(panelMd, XmippWindowUtil.getConstraints(gbc, 0, 1, 2));

		// Checkbox to save images
//		chbImg = new JCheckBox("Save images", false);
//		chbImg.addActionListener(this);
//		group.add(chbImg, XmippWindowUtil.getConstraints(gbc, 0, 2));

		//createImageOptions();
		//group.add(panelImg, XmippWindowUtil.getConstraints(gbc, 0, 3, 2));

		gbc.anchor = GridBagConstraints.PAGE_START;
		panel.add(group, XmippWindowUtil.getConstraints(gbc, 0, 0, 2));
		// Change default Ok text button
		btnOkText = "Save";
		pack();
	}// function initComponents

	protected JPanel createBrowse(BrowseField browse)
	{
		JPanel panel = new JPanel();
		browse.tb = new JTextField(30);
		panel.add(browse.tb);
		browse.btn = XmippWindowUtil.getIconButton("folderopen.gif", this);
		panel.add(browse.btn);
		return panel;
	}

	protected void createMdOptions()
	{
		GridBagConstraints gbc = new GridBagConstraints();
		// Metadata options panel
		panelMd = new JPanel(new GridBagLayout());
		// panelMd.setBackground(Color.blue);
		browseMd = new BrowseField();
		JPanel panelBrowse = createBrowse(browseMd);
		panelMd.add(panelBrowse, XmippWindowUtil.getConstraints(gbc, 0, 0, 3));
		chbDiscard = new JCheckBox("Save Active Metadata Only", false);
		if(!isselection)
			panelMd.add(chbDiscard, XmippWindowUtil.getConstraints(gbc, 0, 1, 2));
		rbMdOverride = new JRadioButton("Overwrite file");
		panelMd.add(rbMdOverride, XmippWindowUtil.getConstraints(gbc, 0, 2));
		gbc.anchor = GridBagConstraints.WEST;
		rbMdAppend = new JRadioButton("Replace/Append", true);
		panelMd.add(rbMdAppend, XmippWindowUtil.getConstraints(gbc, 1, 2));
		ButtonGroup group = new ButtonGroup();
		group.add(rbMdOverride);
		group.add(rbMdAppend);

	}

	public boolean isOverwrite()
	{
		return rbMdOverride.isSelected();
	}

//	protected void createImageOptions()
//	{
//		GridBagConstraints gbc = new GridBagConstraints();
//		// gbc.insets = new Insets(0, 0, 0, 5);
//		// Image options panel
//		panelImg = new JPanel(new GridBagLayout());
//		// panelImg.setBackground(Color.red);
//		panelImg.setVisible(false);
//		gbc.anchor = GridBagConstraints.EAST;
//		panelImg.add(new JLabel("Label  "), XmippWindowUtil.getConstraints(gbc, 0, 0));
//		cbLabel = new JComboBox();
//		for (ColumnInfo ci : data.labels)
//			if (ci.allowRender)
//				cbLabel.addItem(ci.labelName);
//		gbc.anchor = GridBagConstraints.WEST;
//		panelImg.add(cbLabel, XmippWindowUtil.getConstraints(gbc, 1, 0));
//		browseImg = new BrowseField();
//		JPanel panelBrowse = createBrowse(browseImg);
//		panelImg.add(new JLabel("Output stack filename:"), XmippWindowUtil.getConstraints(gbc, 0, 1, 5));
//		panelImg.add(panelBrowse, XmippWindowUtil.getConstraints(gbc, 0, 2, 5));
//		rbStack = new JRadioButton("Stack", true);
//		panelImg.add(rbStack, XmippWindowUtil.getConstraints(gbc, 2, 0));
//		rbIndependent = new JRadioButton("Separate images");
//		// gbc.anchor = GridBagConstraints.WEST;
//		panelImg.add(rbIndependent, XmippWindowUtil.getConstraints(gbc, 3, 0));
//		ButtonGroup group = new ButtonGroup();
//		group.add(rbStack);
//		group.add(rbIndependent);
//	}// function createImageOptions

	/** Reset all controls values to initial values */
	public void setInitialValues()
	{
		browseMd.tb.setText("");
//		browseImg.tb.setText("");
		rbMdAppend.setSelected(true);
//		chbImg.setSelected(false);
//		rbStack.setSelected(true);

	}

	@Override
	public void handleActionPerformed(ActionEvent evt)
	{
		Object obj = evt.getSource();
		// Handle open button action.
		if (obj == browseMd.btn)
			browseFile(browseMd.tb);
		else if (obj == browseImg.btn)
			browseFile(browseImg.tb);
//		else if (obj == chbImg)
//			panelImg.setVisible(chbImg.isSelected());
		else if (obj == chbMd)
			panelMd.setVisible(chbMd.isSelected());

	}

	private void browseFile(JTextField tb)
	{
		int returnVal = fc.showOpenDialog(this);

		if (returnVal == XmippFileChooser.APPROVE_OPTION)
		{
			File file = fc.getSelectedFile();
			String text = (block == null)? file.getPath(): block + "@" + file.getPath(); 
			tb.setText(text);
		}
	}

	public class BrowseField
	{
		public JTextField tb;
		public JButton btn;
	}

	/** Getters and setters */
	public String getMdFilename()
	{
		return browseMd.tb.getText();
	}

	public void setMdFilename(String file)
	{
		if (file != null)
		{
			if (file.contains("@"))
				block = file.substring(0, file.lastIndexOf("@"));
			browseMd.tb.setText(file);
		}
		else
		{
			block = null;
			browseMd.tb.setText("");
		}
			
	}

	public boolean isAppendMode()
	{
		return rbMdAppend.isSelected();
	}

	public void setAppendMode(boolean value)
	{
		rbMdAppend.setSelected(value);
	}

	public boolean doSaveImages()
	{
		return false;
		//return chbImg.isSelected();
	}

	public boolean isOutputIndependent()
	{
		return rbIndependent.isSelected();
	}

	public String getOutput()
	{
		return browseImg.tb.getText();
	}

	public int getImageLabel()
	{
		try
		{
			return MetaData.str2Label(cbLabel.getSelectedItem().toString());
		}
		catch (Exception e)
		{
			e.printStackTrace();
			return MDLabel.MDL_UNDEFINED;
		}
	}

	public boolean saveActiveMetadataOnly()
	{
		return chbDiscard.isSelected();
	}
}// class SaveJDialog
