package particlepicker;

import ij.IJ;
import ij.ImageJ;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;

import particlepicker.tiltpair.gui.TiltPairParticlesJDialog;
import particlepicker.training.gui.TrainingPickerJFrame;
import particlepicker.training.model.TrainingParticle;
import particlepicker.training.model.TrainingPicker;

public abstract class ParticlePickerJFrame extends JFrame
{

	protected ParticlesJDialog particlesdialog;
	private String tool = "Particle Picker Tool";
	protected JMenuItem ijmi;
	
	public ParticlePickerJFrame()
	{
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
	}

	public double getMagnification()
	{
		double scaled = getFamily().getSize()/500.f;
		return 1 - scaled;
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
	
	public abstract boolean isShapeSelected(Shape shape);
	
	public Tool getTool()
	{

		if (IJ.getInstance() == null)
			return Tool.PICKER;
		if (IJ.getToolName().equalsIgnoreCase(tool))
			return Tool.PICKER;
		return Tool.IMAGEJ;
	}
	

}
