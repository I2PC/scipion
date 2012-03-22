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
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import xmipp.utils.XmippWindowUtil;
import xmipp.utils.XmippMessage;

public class ImportParticlesFromProjectJDialog extends JDialog
{

	ParticlePickerJFrame parent;
	private JRadioButton xmipp24rb;
	private JRadioButton xmipp30rb;
	private JRadioButton emanrb;
	private JPanel sourcepn;
	private JTextField sourcetf;
	private ButtonGroup formatgroup;
	public Format format;

	public ImportParticlesFromProjectJDialog(ParticlePickerJFrame parent, boolean modal)
	{
		super(parent, modal);
		setResizable(false);
		this.parent = parent;
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Import from Project...");
		setLayout(new GridBagLayout());
		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(5, 5, 5, 5);
		constraints.anchor = GridBagConstraints.WEST;
		initSourcePane();
		add(new JLabel("Format:"), XmippWindowUtil.getConstraints(constraints, 0, 0, 1));
		add(sourcepn, XmippWindowUtil.getConstraints(constraints, 1, 0, 2));
		add(new JLabel("Source:"), XmippWindowUtil.getConstraints(constraints, 0, 1, 1));
		sourcetf = new JTextField(20);
		add(sourcetf, XmippWindowUtil.getConstraints(constraints, 1, 1, 1));
		JButton browsebt = new JButton("Browse");
		add(browsebt, XmippWindowUtil.getConstraints(constraints, 2, 1, 1));
		browsebt.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				browseDirectory();
			}

		});
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
					JOptionPane.showMessageDialog(ImportParticlesFromProjectJDialog.this, ex.getMessage());
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
		add(actionspn, XmippWindowUtil.getConstraints(constraints, 0, 2, 3));
		pack();
		XmippWindowUtil.setLocation(0.8f, 0.5f, this);
		setVisible(true);

	}

	private void browseDirectory()
	{
		XmippFileChooser fc = new XmippFileChooser(new File(parent.getParticlePicker().getOutputDir()));
		fc.setFileSelectionMode(XmippFileChooser.DIRECTORIES_ONLY);
		int returnVal = fc.showOpenDialog(ImportParticlesFromProjectJDialog.this);

		try
		{
			if (returnVal == XmippFileChooser.APPROVE_OPTION)
			{
				String path = fc.getSelectedFile().getAbsolutePath();
				sourcetf.setText(path);
			}
		}
		catch (Exception ex)
		{
			JOptionPane.showMessageDialog(ImportParticlesFromProjectJDialog.this, ex.getMessage());
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

		emanrb = new JRadioButton(Format.Eman.toString());
		emanrb.addItemListener(formatlistener);

		sourcepn.add(xmipp24rb);
		sourcepn.add(xmipp30rb);
		sourcepn.add(emanrb);

		formatgroup.add(xmipp24rb);
		formatgroup.add(xmipp30rb);
		formatgroup.add(emanrb);

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

		String projectdir = sourcetf.getText();
		if (projectdir == null || projectdir.equals(""))
			throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("directory"));
		switch (format)
		{
		case Xmipp24:
			parent.importParticlesXmipp24(projectdir);
			break;
		case Xmipp30:
			parent.importParticlesXmipp30(projectdir);
			break;
		case Eman:
			parent.importParticlesEman(projectdir);
			break;

		}

	}
}
