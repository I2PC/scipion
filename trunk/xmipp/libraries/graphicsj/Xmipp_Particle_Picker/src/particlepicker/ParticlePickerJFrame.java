package particlepicker;

import ij.IJ;
import ij.ImageJ;

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
import java.text.NumberFormat;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
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

import particlepicker.tiltpair.gui.TiltPairParticlesJDialog;
import particlepicker.training.gui.TrainingPickerJFrame;
import particlepicker.training.model.TrainingParticle;
import particlepicker.training.model.TrainingPicker;

public abstract class ParticlePickerJFrame extends JFrame implements ActionListener
{

	protected ParticlesJDialog particlesdialog;
	private String tool = "Particle Picker Tool";
	protected JMenuItem ijmi;
	protected JCheckBox circlechb;
	protected JCheckBox rectanglechb;
	protected JFormattedTextField sizetf;
	protected JCheckBox centerchb;
	protected JPanel symbolpn;
	protected JMenuItem savemi;
	protected JMenuItem hcontentsmi;
	protected JMenuItem pmi;
	protected JMenu filtersmn;
	protected String activemacro;
	protected JSlider sizesl;
	protected JPanel sizepn;
	
	
	public ParticlePickerJFrame()
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
		
		
		
		ijmi = new JMenuItem("ImageJ");
		ijmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				if (IJ.getInstance() == null)
				{

					new ImageJ();
					IJ.run("Install...", "install=" + TrainingPicker.getXmippPath("external/imagej/macros/ParticlePicker.txt"));
					IJ.setTool(tool);
				}
				// IJ.getInstance().setVisible(true);
			}
		});
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
		hcontentsmi = new JMenuItem("Help Contents...");
		hcontentsmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				try
				{
					WindowUtils.openURI("http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParticlePicker");
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
		filtersmn = new JMenu("Filters");
		JCheckBoxMenuItem fftbpf = new JCheckBoxMenuItem("Bandpass Filter...");
		filtersmn.add(fftbpf);
		fftbpf.addActionListener(this);
		JCheckBoxMenuItem admi = new JCheckBoxMenuItem("Anisotropic Diffusion...");
		filtersmn.add(admi);
		admi.addActionListener(new ActionListener()
		{
			
			@Override
			public void actionPerformed(ActionEvent e)
			{
				activemacro = "8-bit";
				IJ.run(activemacro);
				activemacro = ((JCheckBoxMenuItem) e.getSource()).getText();
				IJ.run(activemacro);
			}
		});
		JCheckBoxMenuItem msmi = new JCheckBoxMenuItem("Mean Shift");
		filtersmn.add(msmi);
		msmi.addActionListener(this);
		JCheckBoxMenuItem sbmi = new JCheckBoxMenuItem("Subtract Background...");
		filtersmn.add(sbmi);
		sbmi.addActionListener(this);
		JCheckBoxMenuItem gbmi = new JCheckBoxMenuItem("Gaussian Blur...");
		filtersmn.add(gbmi);
		gbmi.addActionListener(this);
		JCheckBoxMenuItem bcmi = new JCheckBoxMenuItem("Brightness/Contrast...");
		filtersmn.add(bcmi);
		bcmi.addActionListener(this);
		
	}
	
	@Override
	public void actionPerformed(ActionEvent e)
	{
		try
		{
			activemacro = ((JCheckBoxMenuItem) e.getSource()).getText();
			IJ.run(activemacro);
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
				if(ParticlePickerJFrame.this instanceof TrainingPickerJFrame)
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
			if(particlesdialog != null)
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
	
	
	public Tool getTool()
	{

		if (IJ.getInstance() == null)
			return Tool.PICKER;
		if (IJ.getToolName().equalsIgnoreCase(tool))
			return Tool.PICKER;
		return Tool.IMAGEJ;
	}
	
	public boolean isPickingAvailable(MouseEvent e)
	{
		if(getTool() != Tool.PICKER)
			return false;
		if(SwingUtilities.isRightMouseButton(e))
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
			getCanvas().repaint();
		}
	}
	
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

}
