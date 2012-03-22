package xmipp.particlepicker;

import ij.CommandListener;
import ij.Executer;
import ij.IJ;
import ij.ImageJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.plugin.frame.Recorder;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import xmipp.utils.XmippFileChooser;
import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;

import xmipp.ij.commons.Tool;
import xmipp.jni.Filename;
import xmipp.particlepicker.tiltpair.gui.TiltPairParticlesJDialog;
import xmipp.particlepicker.training.gui.TrainingPickerJFrame;
import xmipp.particlepicker.training.model.TrainingParticle;
import xmipp.particlepicker.training.model.TrainingPicker;
import xmipp.utils.XmippWindowUtil;
import xmipp.utils.XmippMessage;

public abstract class ParticlePickerJFrame extends JFrame implements ActionListener
{

	protected ParticlesJDialog particlesdialog;
	
	protected JMenuItem ijmi;
	protected JCheckBox circlechb;
	protected JCheckBox rectanglechb;
	protected JFormattedTextField sizetf;
	protected JCheckBox centerchb;
	protected JPanel symbolpn;
	protected JMenuItem savemi;
	protected JMenuItem hcontentsmi;
	protected JMenuItem pmi;
	protected JMenuItem importfpmi;
	protected JMenuItem exportmi;
	protected JMenu filtersmn;
	protected String activefilter;
	protected JSlider sizesl;
	protected JPanel sizepn;
	private String command;
	private List<JCheckBoxMenuItem> mifilters;
	protected JMenu filemn;
	protected JMenuItem importffmi;

	public ParticlePickerJFrame(ParticlePicker picker)
	{

		addWindowListener(new WindowAdapter()
		{
			public void windowClosing(WindowEvent winEvt)
			{
				if (getParticlePicker().isChanged())
				{
					int result = JOptionPane.showConfirmDialog(ParticlePickerJFrame.this, "Save changes before closing?", "Message", JOptionPane.YES_NO_OPTION);
					if (result == JOptionPane.OK_OPTION)
						saveChanges();
				}
				System.exit(0);
			}
		});

		Recorder.record = true;

		// detecting if a command is thrown by ImageJ
		Executer.addCommandListener(new CommandListener()
		{
			public String commandExecuting(String command)
			{
				ParticlePickerJFrame.this.command = command;
				return command;

			}
		});
		ImagePlus.addImageListener(new ImageListener()
		{

			@Override
			public void imageUpdated(ImagePlus arg0)
			{
				if (command != null)
				{
					String options = "";
					if (Recorder.getCommandOptions() != null)
						options = Recorder.getCommandOptions();
					if (!getParticlePicker().isFilterSelected(command))
						getParticlePicker().addFilter(command, options);
					command = null;
				}
			}

			@Override
			public void imageOpened(ImagePlus arg0)
			{
				// TODO Auto-generated method stub

			}

			@Override
			public void imageClosed(ImagePlus arg0)
			{
				// TODO Auto-generated method stub

			}
		});
		filemn = new JMenu("File");
		savemi = new JMenuItem("Save");
		savemi.setMnemonic('S');
		savemi.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.CTRL_DOWN_MASK));
		savemi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				saveChanges();
				JOptionPane.showMessageDialog(ParticlePickerJFrame.this, "Data saved successfully");
				((JMenuItem) e.getSource()).setEnabled(false);
			}
		});
		filemn.add(savemi);
		importfpmi = new JMenuItem("Import from Project...");
		filemn.add(importfpmi);
		importfpmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				new ImportParticlesFromProjectJDialog(ParticlePickerJFrame.this, true);
			}
		});
		
		importffmi = new JMenuItem("Import from Files...");
		filemn.add(importffmi);
		importffmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				displayImportDialog();
				
			}
		});
		
		exportmi = new JMenuItem("Export Particles...");
		filemn.add(exportmi);
		exportmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				XmippFileChooser fc = new XmippFileChooser();
				int returnVal = fc.showOpenDialog(ParticlePickerJFrame.this);

				try
				{
					if (returnVal == XmippFileChooser.APPROVE_OPTION)
					{
						File file = fc.getSelectedFile();
						getParticlePicker().exportParticles(getFamily(), file.getAbsolutePath());
						JOptionPane.showMessageDialog(ParticlePickerJFrame.this, "Export successful");
					}
				}
				catch (Exception ex)
				{
					JOptionPane.showMessageDialog(ParticlePickerJFrame.this, ex.getMessage());
				}
			}
		});

		ijmi = new JMenuItem("ImageJ");
		ijmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				if (IJ.getInstance() == null)
				{

					new ImageJ();
					IJ.run("Install...", "install=" + Filename.getXmippPath("java/src/xmipp/ij/XmippMacros.txt"));
					IJ.setTool(Tool.getTool(Tool.PICKER));
				}
			}
		});
		
		hcontentsmi = new JMenuItem("Help Contents...");
		hcontentsmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				try
				{
					XmippWindowUtil.openURI("http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParticlePicker");
				}
				catch (Exception ex)
				{
					JOptionPane.showMessageDialog(ParticlePickerJFrame.this, ex.getMessage());
				}
			}
		});
		pmi = new JMenuItem("Particles");
		pmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				loadParticles();
			}
		});

		mifilters = new ArrayList<JCheckBoxMenuItem>();
		filtersmn = new JMenu("Filters");
		filtersmn.addMenuListener(new MenuListener()
		{

			@Override
			public void menuCanceled(MenuEvent arg0)
			{
				// TODO Auto-generated method stub

			}

			@Override
			public void menuDeselected(MenuEvent arg0)
			{
				// TODO Auto-generated method stub

			}

			@Override
			public void menuSelected(MenuEvent arg0)
			{
				for (JCheckBoxMenuItem mi : mifilters)
					mi.setSelected(getParticlePicker().isFilterSelected(mi.getText()));

			}
		});

		addFilterMenuItem("Bandpass Filter...", true, picker);
		JCheckBoxMenuItem admi = addFilterMenuItem("Anisotropic Diffusion...", false, picker);
		admi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				activefilter = "8-bit";
				IJ.run(activefilter);
				activefilter = ((JCheckBoxMenuItem) e.getSource()).getText();
				IJ.run(activefilter);
			}
		});
		addFilterMenuItem("Mean Shift", true, picker);
		addFilterMenuItem("Subtract Background...", true, picker);
		addFilterMenuItem("Gaussian Blur...", true, picker);
		addFilterMenuItem("Brightness/Contrast...", true, picker);
	}

	protected abstract void displayImportDialog();

	private JCheckBoxMenuItem addFilterMenuItem(String command, boolean defaultlistener, ParticlePicker picker)
	{
		JCheckBoxMenuItem mi = new JCheckBoxMenuItem(command);
		mifilters.add(mi);
		mi.setSelected(picker.isFilterSelected(command));
		if (defaultlistener)
			mi.addActionListener(this);
		filtersmn.add(mi);
		return mi;
	}

	@Override
	public void actionPerformed(ActionEvent e)
	{
		try
		{
			JCheckBoxMenuItem item = (JCheckBoxMenuItem) e.getSource();
			activefilter = item.getText();
			if (item.isSelected())// filter added, will be registered by picker
									// with options if needed
				IJ.run(activefilter);
			else
				// filter removed
				getParticlePicker().removeFilter(activefilter);
			setChanged(true);
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
			JOptionPane.showMessageDialog(this, ex.getMessage());
		}

	}

	protected abstract void saveChanges();

	public int getSide(int size)
	{
		return 100;
	}

	public abstract Family getFamily();

	public abstract ParticlePickerCanvas getCanvas();

	public void loadParticles()
	{
		try
		{
			if (particlesdialog == null)
				if (ParticlePickerJFrame.this instanceof TrainingPickerJFrame)
					particlesdialog = new ParticlesJDialog(ParticlePickerJFrame.this);
				else
					particlesdialog = new TiltPairParticlesJDialog(ParticlePickerJFrame.this);
			else
			{

				particlesdialog.loadParticles(true);
				particlesdialog.setVisible(true);
			}
		}
		catch (Exception ex)
		{
			JOptionPane.showMessageDialog(ParticlePickerJFrame.this, ex.getMessage());
			if (particlesdialog != null)
				particlesdialog.close();
			particlesdialog = null;
		}
	}

	public void updateMicrographsModel()
	{
		if (particlesdialog != null)
			loadParticles();
	}

	public ParticlesJDialog getParticlesJDialog()
	{
		return particlesdialog;
	}

	public abstract Micrograph getMicrograph();

	public abstract List<? extends TrainingParticle> getParticles();



	public boolean isPickingAvailable(MouseEvent e)
	{
		if (getCanvas().getTool() != Tool.PICKER)
			return false;
		if (SwingUtilities.isRightMouseButton(e))
			return false;
		return true;
	}

	protected void initSymbolPane()
	{

		symbolpn = new JPanel(new FlowLayout(FlowLayout.LEFT));
		symbolpn.add(new JLabel("Symbol:"));
		// symbolpn.setBorder(BorderFactory.createTitledBorder("Symbol"));
		ShapeItemListener shapelistener = new ShapeItemListener();

		circlechb = new JCheckBox(Shape.Circle.toString());
		circlechb.setSelected(true);
		circlechb.addItemListener(shapelistener);

		rectanglechb = new JCheckBox(Shape.Rectangle.toString());
		rectanglechb.setSelected(true);
		rectanglechb.addItemListener(shapelistener);

		centerchb = new JCheckBox(Shape.Center.toString());
		centerchb.setSelected(true);
		centerchb.addItemListener(shapelistener);

		symbolpn.add(circlechb);
		symbolpn.add(rectanglechb);
		symbolpn.add(centerchb);
	}

	class ShapeItemListener implements ItemListener
	{
		@Override
		public void itemStateChanged(ItemEvent e)
		{
			changeShapes();
		}
	}

	public abstract void changeShapes();

	public boolean isShapeSelected(Shape shape)
	{
		switch (shape)
		{
		case Rectangle:
			return rectanglechb.isSelected();
		case Circle:
			return circlechb.isSelected();
		case Center:
			return centerchb.isSelected();
			// case OnlyLast:
			// return onlylastchb.isSelected();
		}
		return false;
	}

	public abstract ParticlePicker getParticlePicker();

	public abstract void setChanged(boolean changed);

	protected void initSizePane()
	{
		sizepn = new JPanel();

		int size = getFamily().getSize();
		sizepn.add(new JLabel("Size:"));
		sizesl = new JSlider(0, 1000, size);
		sizesl.setPaintTicks(true);
		sizesl.setMajorTickSpacing(100);
		int height = (int) sizesl.getPreferredSize().getHeight();

		sizesl.setPreferredSize(new Dimension(100, height));
		sizepn.add(sizesl);

		sizetf = new JFormattedTextField(NumberFormat.getIntegerInstance());
		sizetf.setColumns(3);
		sizetf.setText(Integer.toString(size));
		sizepn.add(sizetf);
		sizetf.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				int size = ((Number) sizetf.getValue()).intValue();
				switchSize(size);
			}
		});

		sizesl.addChangeListener(new ChangeListener()
		{

			@Override
			public void stateChanged(ChangeEvent e)
			{
				int size = sizesl.getValue();
				switchSize(size);
			}
		});

	}

	public void switchSize(int size)
	{
		sizetf.setText(Integer.toString(size));
		sizesl.setValue(size);
		getCanvas().repaint();
		getFamily().setSize(size);
		if (particlesdialog != null)
		{
			for (TrainingParticle p : getParticles())
				p.resetParticleCanvas();
			loadParticles();
		}
		setChanged(true);
	}
	
	protected void importParticlesXmipp30(String path)
	{
		getParticlePicker().importParticlesXmipp30(getFamily(), path);
		setChanged(true);
		getCanvas().repaint();
		updateMicrographsModel();
	}
	
	
	public  void importParticlesXmipp24(String projectdir)
	{
		getParticlePicker().importParticlesFromXmipp24Project(getFamily(), projectdir);
		setChanged(true);
		getCanvas().repaint();
		updateMicrographsModel();
		
	}

	public  void importParticlesEman(String path)
	{
		throw new UnsupportedOperationException(XmippMessage.getNotImplementedYetMsg());
		
	}



	


}
