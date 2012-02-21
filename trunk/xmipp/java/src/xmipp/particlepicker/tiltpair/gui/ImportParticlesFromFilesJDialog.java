package xmipp.particlepicker.tiltpair.gui;

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

public class ImportParticlesFromFilesJDialog extends JDialog
{

	private TiltPairPickerJFrame parent;
	private JRadioButton xmipp24rb;
	private JRadioButton xmipp30rb;
	private JPanel sourcepn;
	private JTextField untiltedtf;
	private ButtonGroup formatgroup;
	public Format format;
	private JTextField tiltedtf;
	private JButton browsetbt;
	private JButton browseubt;
	private JFileChooser fc;

	public ImportParticlesFromFilesJDialog(TiltPairPickerJFrame parent, boolean modal)
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
		add(new JLabel("Untilted:"), WindowUtil.getConstraints(constraints, 0, 1, 1));
		untiltedtf = new JTextField(20);

		add(untiltedtf, WindowUtil.getConstraints(constraints, 1, 1, 1));
		BrowseListener bl = new BrowseListener();
		browseubt = new JButton("Browse");
		add(browseubt, WindowUtil.getConstraints(constraints, 2, 1, 1));
		browseubt.addActionListener(bl);
		add(new JLabel("Tilted:"), WindowUtil.getConstraints(constraints, 0, 2, 1));
		tiltedtf = new JTextField(20);
		add(tiltedtf, WindowUtil.getConstraints(constraints, 1, 2, 1));
		browsetbt = new JButton("Browse");
		add(browsetbt, WindowUtil.getConstraints(constraints, 2, 2, 1));
		browsetbt.addActionListener(bl);

		fc = new JFileChooser();
		setFile(untiltedtf, parent.getMicrograph().getPosFileFromXmipp24());
		setFile(tiltedtf, parent.getMicrograph().getTiltedMicrograph().getPosFileFromXmipp24());

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
					JOptionPane.showMessageDialog(ImportParticlesFromFilesJDialog.this, ex.getMessage());
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
		if(filepath == null)
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
			JTextField tf = (e.getSource().equals(browseubt)) ? untiltedtf : tiltedtf;
			browseDirectory(tf);
		}

	}

	private void browseDirectory(JTextField tf)
	{

		int returnVal = fc.showOpenDialog(ImportParticlesFromFilesJDialog.this);

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
			JOptionPane.showMessageDialog(ImportParticlesFromFilesJDialog.this, ex.getMessage());
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

		String ufile = untiltedtf.getText();
		String tfile = tiltedtf.getText();
		if (ufile == null || ufile.equals("") || ufile == null || ufile.equals(""))
			throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("files"));
		switch (format)
		{
		case Xmipp24:
			parent.importParticlesFromXmipp24Files(ufile, tfile);
			break;
		case Xmipp30:
			parent.importParticlesFromXmipp30Files(ufile, tfile);
			break;

		}

	}
}
