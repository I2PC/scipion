package xmipp.particlepicker.training.gui;

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
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;

import xmipp.particlepicker.Format;
import xmipp.particlepicker.ParticlePickerJFrame;
import xmipp.utils.XmippWindowUtil;
import xmipp.utils.XmippMessage;

public class ImportParticlesFromFileJDialog extends JDialog {

	private TrainingPickerJFrame parent;
	private JRadioButton xmipp24rb;
	private JRadioButton xmipp30rb;
	private JPanel sourcepn;
	private JTextField filetf;
	private ButtonGroup formatgroup;
	public Format format;
	private JTextField tiltedtf;
	private JButton browsetbt;
	private JButton browseubt;
	private XmippFileChooser fc;
	private JRadioButton xmippemanrb;

	public ImportParticlesFromFileJDialog(TrainingPickerJFrame parent,
			boolean modal) {
		super(parent, modal);
		setResizable(false);
		this.parent = parent;
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Import from Files...");
		setLayout(new GridBagLayout());
		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(5, 5, 5, 5);
		constraints.anchor = GridBagConstraints.WEST;
		initSourcePane();
		add(new JLabel("Format:"),
				XmippWindowUtil.getConstraints(constraints, 0, 0, 1));
		add(sourcepn, XmippWindowUtil.getConstraints(constraints, 1, 0, 2));
		add(new JLabel("File:"),
				XmippWindowUtil.getConstraints(constraints, 0, 1, 1));
		filetf = new JTextField(20);

		add(filetf, XmippWindowUtil.getConstraints(constraints, 1, 1, 1));
		BrowseListener bl = new BrowseListener();
		browseubt = XmippWindowUtil.getIconButton("folderopen.gif", null);
		add(browseubt, XmippWindowUtil.getConstraints(constraints, 2, 1, 1));
		browseubt.addActionListener(bl);

		fc = new XmippFileChooser(new File(parent.getParticlePicker().getOutputDir()));
		setFile(filetf, parent.getMicrograph().getPosFileFromXmipp24());

		JPanel actionspn = new JPanel(new FlowLayout(FlowLayout.RIGHT));
		JButton okbt = XmippWindowUtil.getTextButton("OK",
				new ActionListener() {
					@Override
					public void actionPerformed(ActionEvent arg0) {
						try {
							importParticles();
							setVisible(false);
							dispose();
						} catch (Exception ex) {
							JOptionPane.showMessageDialog(
									ImportParticlesFromFileJDialog.this,
									ex.getMessage());
						}
					}
				});
		JButton cancelbt = XmippWindowUtil.getTextButton("Cancel",
				new ActionListener() {
					@Override
					public void actionPerformed(ActionEvent arg0) {
						setVisible(false);
						dispose();
					}
				});
		actionspn.add(okbt);
		actionspn.add(cancelbt);
		add(actionspn, XmippWindowUtil.getConstraints(constraints, 0, 3, 3));

		pack();
		XmippWindowUtil.setLocation(0.8f, 0.5f, this);
		setVisible(true);

	}

	private void setFile(JTextField tf, String filepath) {
		if (filepath == null)
			return;
		File file = new File(filepath);
		if (file.exists()) {
			tf.setText(filepath);
			fc.setSelectedFile(file);
		}
	}

	class BrowseListener implements ActionListener {

		@Override
		public void actionPerformed(ActionEvent e) {
			JTextField tf = (e.getSource().equals(browseubt)) ? filetf
					: tiltedtf;
			browseDirectory(tf);
		}

	}

	private void browseDirectory(JTextField tf) {

		int returnVal = fc.showOpenDialog(ImportParticlesFromFileJDialog.this);

		try {
			if (returnVal == XmippFileChooser.APPROVE_OPTION) {
				String path = fc.getSelectedFile().getAbsolutePath();
				tf.setText(path);
			}
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(ImportParticlesFromFileJDialog.this,
					ex.getMessage());
		}

	}

	protected void initSourcePane() {

		sourcepn = new JPanel(new FlowLayout(FlowLayout.LEFT));

		FormatItemListener formatlistener = new FormatItemListener();

		formatgroup = new ButtonGroup();
		
		xmipp24rb = new JRadioButton(Format.Xmipp24.toString());
		xmipp24rb.setSelected(true);
		format = Format.Xmipp24;
		xmipp24rb.addItemListener(formatlistener);
		formatgroup.add(xmipp24rb);
		
		xmipp30rb = new JRadioButton(Format.Xmipp30.toString());
		xmipp30rb.addItemListener(formatlistener);
		formatgroup.add(xmipp30rb);
		
		xmippemanrb = new JRadioButton(Format.Eman.toString());
		xmippemanrb.addItemListener(formatlistener);
		formatgroup.add(xmippemanrb);
		
		sourcepn.add(xmipp24rb);
		sourcepn.add(xmipp30rb);
		sourcepn.add(xmippemanrb);

	}

	class FormatItemListener implements ItemListener {

		@Override
		public void itemStateChanged(ItemEvent e) {
			JRadioButton formatrb = (JRadioButton) e.getSource();
			format = Format.valueOf(formatrb.getText());
		}
	}

	private void importParticles() {

		String file = filetf.getText();
		if (file == null || file.equals(""))
			throw new IllegalArgumentException(
					XmippMessage.getEmptyFieldMsg("file"));
		parent.importParticlesFromFile(format, file);

	}
}
