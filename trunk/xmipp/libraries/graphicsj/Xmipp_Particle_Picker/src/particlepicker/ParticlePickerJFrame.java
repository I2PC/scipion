package particlepicker;

import java.util.List;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

import particlepicker.training.gui.TrainingPickerJFrame;
import particlepicker.training.model.Family;
import particlepicker.training.model.Micrograph;
import particlepicker.training.model.TrainingParticle;

public abstract class ParticlePickerJFrame extends JFrame
{

	protected ParticlesJDialog particlesdialog;

	public double getMagnification()
	{
		if (getFamily().getSize() > 100)
			return 0.9;
		return 1;
	}

	public abstract Family getFamily();

	public abstract ParticlePickerCanvas getCanvas();

	public void loadParticles()
	{
		try
		{
			if (particlesdialog == null)
				particlesdialog = new ParticlesJDialog(ParticlePickerJFrame.this);
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

}
