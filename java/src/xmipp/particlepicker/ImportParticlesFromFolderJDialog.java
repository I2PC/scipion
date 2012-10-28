package xmipp.particlepicker;

import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JDialog;
import xmipp.utils.XmippFileChooser;

import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;

import xmipp.utils.XmippDialog;
import xmipp.utils.XmippWindowUtil;
import xmipp.utils.XmippMessage;

public class ImportParticlesFromFolderJDialog extends XmippDialog {

	ParticlePickerJFrame parent;
	private JRadioButton xmipp24rb;
	private JRadioButton xmipp30rb;
	private JRadioButton emanrb;
	private JPanel sourcepn;
	private JTextField sourcetf;
	private JButton browsebt;
	private JComboBox jcbFormat;
	private ButtonGroup formatgroup;
	public Format format;
	
	private static String[] FormatStrings = {"Automatic", "Xmipp 2.4", "Xmipp 3.0", "Eman"};
	private static Format[] FormatList = {Format.Auto, Format.Xmipp24, Format.Xmipp30, Format.Eman};

	public ImportParticlesFromFolderJDialog(JFrame parent) {
		super(parent, "Import from Project...", true);
		this.parent = (ParticlePickerJFrame) parent;
	}
	
	@Override
	protected void createContent(JPanel panel) {
		setResizable(false);
		panel.setLayout(new GridBagLayout());
		GridBagConstraints gbc = new GridBagConstraints();
		gbc.insets = new Insets(5, 5, 5, 5);
		gbc.anchor = GridBagConstraints.WEST;
		panel.add(new JLabel("Format:"),
				XmippWindowUtil.getConstraints(gbc, 0, 0, 1));
		panel.add(sourcepn, XmippWindowUtil.getConstraints(gbc, 1, 0, 2));
		panel.add(new JLabel("Source:"),
				XmippWindowUtil.getConstraints(gbc, 0, 1, 1));
		sourcetf = new JTextField(20);
		panel.add(sourcetf, XmippWindowUtil.getConstraints(gbc, 1, 1, 1));
		browsebt = XmippWindowUtil.getIconButton("folderopen.gif", this);
		
		/** Create a combobox with possible formats */
		jcbFormat = new JComboBox(FormatStrings);
     	jcbFormat.addActionListener(this);
//				new ActionListener() {
//					@Override
//					public void actionPerformed(ActionEvent e) {
//						browseDirectory();
//					}
//				});
		panel.add(browsebt, XmippWindowUtil.getConstraints(gbc, 2, 1, 1));
//		JPanel actionspn = new JPanel(new FlowLayout(FlowLayout.RIGHT));
//		JButton okbt = XmippWindowUtil.getTextButton("OK",
//				new ActionListener() {
//					@Override
//					public void actionPerformed(ActionEvent arg0) {
//						try {
//							importParticles();
//							setVisible(false);
//							dispose();
//						} catch (Exception ex) {
//							JOptionPane.showMessageDialog(
//									ImportParticlesFromFolderJDialog.this,
//									ex.getMessage());
//						}
//					}
//				});
//		JButton cancelbt = XmippWindowUtil.getTextButton("Cancel",
//				new ActionListener() {
//					@Override
//					public void actionPerformed(ActionEvent arg0) {
//						setVisible(false);
//						dispose();
//					}
//				});
//		actionspn.add(okbt);
//		actionspn.add(cancelbt);
//		add(actionspn, XmippWindowUtil.getConstraints(gbc, 0, 2, 3));
//		pack();
//		XmippWindowUtil.setLocation(0.8f, 0.5f, this);
//		setVisible(true);

	}
	
	@Override
	public void handleActionPerformed(ActionEvent evt) {
		Object o = evt.getSource();
		if (o == browsebt)
			browseDirectory();
//		else if (o == jcbLabel)
//			label = labels[jcbLabel.getSelectedIndex()];
//		else if (o == btnReplace)
//			replace(selected_index);
//		else if (o == btnReplaceFind) {
//			replace(selected_index);
//			find();
//		} else if (o == btnReplaceAll)
//			replaceAll();

	}// function actionPerformed	
	

	private void browseDirectory() {
		XmippFileChooser fc = new XmippFileChooser(new File(parent
				.getParticlePicker().getOutputDir()));
		fc.setFileSelectionMode(XmippFileChooser.DIRECTORIES_ONLY);
		int returnVal = fc
				.showOpenDialog(ImportParticlesFromFolderJDialog.this);

		try {
			if (returnVal == XmippFileChooser.APPROVE_OPTION) {
				String path = fc.getSelectedFile().getAbsolutePath();
				sourcetf.setText(path);
			}
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(
					ImportParticlesFromFolderJDialog.this, ex.getMessage());
		}

	}

//	protected void initSourcePane() {
//
//		sourcepn = new JPanel(new FlowLayout(FlowLayout.LEFT));
//
//		FormatItemListener formatlistener = new FormatItemListener();
//
//		formatgroup = new ButtonGroup();
//		xmipp24rb = new JRadioButton(Format.Xmipp24.toString());
//		xmipp24rb.setSelected(true);
//		format = Format.Xmipp24;
//		xmipp24rb.addItemListener(formatlistener);
//
//		xmipp30rb = new JRadioButton(Format.Xmipp30.toString());
//		xmipp30rb.addItemListener(formatlistener);
//
//		emanrb = new JRadioButton(Format.Eman.toString());
//		emanrb.addItemListener(formatlistener);
//
//		sourcepn.add(xmipp24rb);
//		sourcepn.add(xmipp30rb);
//		sourcepn.add(emanrb);
//
//		formatgroup.add(xmipp24rb);
//		formatgroup.add(xmipp30rb);
//		formatgroup.add(emanrb);
//
//	}

//	class FormatItemListener implements ItemListener {
//
//		@Override
//		public void itemStateChanged(ItemEvent e) {
//			JRadioButton formatrb = (JRadioButton) e.getSource();
//			format = Format.valueOf(formatrb.getText());
//		}
//	}

	private void importParticles() {

		
		String dir = sourcetf.getText();
		if (dir == null || dir.equals(""))
			throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("directory"));
		parent.importParticlesFromFolder(format, dir);
			
		

	}
}
