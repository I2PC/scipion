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
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;

import xmipp.particlepicker.Format;
import xmipp.particlepicker.ParticlePickerJFrame;
import xmipp.utils.WindowUtil;
import xmipp.utils.XmippMessage;

public class ImportParticlesFromFileJDialog extends JDialog
{

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
	private JFileChooser fc;

	public ImportParticlesFromFileJDialog(TrainingPickerJFrame parent, boolean modal)
	{
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
		add(new JLabel("Format:"), WindowUtil.getConstraints(constraints, 0, 0, 1));
		add(sourcepn, WindowUtil.getConstraints(constraints, 1, 0, 2));
		add(new JLabel("File:"), WindowUtil.getConstraints(constraints, 0, 1, 1));
		filetf = new JTextField(20);

		add(filetf, WindowUtil.getConstraints(constraints, 1, 1, 1));
		BrowseListener bl = new BrowseListener();
		browseubt = new JButton("Browse");
		add(browseubt, WindowUtil.getConstraints(constraints, 2, 1, 1));
		browseubt.addActionListener(bl);

		fc = new JFileChooser(new File(parent.getParticlePicker().getOutputDir()));
		setFile(filetf, parent.getMicrograph().getPosFileFromXmipp24());

		JPanel actionspn = new JPanel(new FlowLayout(FlowLayout.RIGHT));
		JButton okbt = new JButton("OK");
		okbt.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				try
				{
					importParticles();
					setVisible(false);
					dispose();
				}
				catch (Exception ex)
				{
					JOptionPane.showMessageDialog(ImportParticlesFromFileJDialog.this, ex.getMessage());
				}
			}
		});
		JButton cancelbt = new JButton("Cancel");
		cancelbt.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				setVisible(false);
				dispose();

			}
		});
		actionspn.add(okbt);
		actionspn.add(cancelbt);
		add(actionspn, WindowUtil.getConstraints(constraints, 0, 3, 3));

		pack();
		WindowUtil.setLocation(0.8f, 0.5f, this);
		setVisible(true);

	}

	private void setFile(JTextField tf, String filepath)
	{
		if (filepath == null)
			return;
		File file = new File(filepath);
		if (file.exists())
		{
			tf.setText(filepath);
			fc.setSelectedFile(file);
		}
	}

	class BrowseListener implements ActionListener
	{

		@Override
		public void actionPerformed(ActionEvent e)
		{
			JTextField tf = (e.getSource().equals(browseubt)) ? filetf : tiltedtf;
			browseDirectory(tf);
		}

	}

	private void browseDirectory(JTextField tf)
	{

		int returnVal = fc.showOpenDialog(ImportParticlesFromFileJDialog.this);

		try
		{
			if (returnVal == JFileChooser.APPROVE_OPTION)
			{
				String path = fc.getSelectedFile().getAbsolutePath();
				tf.setText(path);
			}
		}
		catch (Exception ex)
		{
			JOptionPane.showMessageDialog(ImportParticlesFromFileJDialog.this, ex.getMessage());
		}

	}

	protected void initSourcePane()
	{

		sourcepn = new JPanel(new FlowLayout(FlowLayout.LEFT));

		FormatItemListener formatlistener = new FormatItemListener();

		formatgroup = new ButtonGroup();
		xmipp24rb = new JRadioButton(Format.Xmipp24.toString());
		xmipp24rb.setSelected(true);
		format = Format.Xmipp24;
		xmipp24rb.addItemListener(formatlistener);

		xmipp30rb = new JRadioButton(Format.Xmipp30.toString());
		xmipp30rb.addItemListener(formatlistener);

		sourcepn.add(xmipp24rb);
		sourcepn.add(xmipp30rb);

		formatgroup.add(xmipp24rb);
		formatgroup.add(xmipp30rb);

	}

	class FormatItemListener implements ItemListener
	{

		@Override
		public void itemStateChanged(ItemEvent e)
		{
			JRadioButton formatrb = (JRadioButton) e.getSource();
			format = Format.valueOf(formatrb.getText());
		}
	}

	private void importParticles()
	{

		String file = filetf.getText();
		if (file == null || file.equals(""))
			throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("file"));
		switch (format)
		{
		case Xmipp30:
			parent.importParticlesFromXmipp30File(file);
			break;

		}

	}
}
