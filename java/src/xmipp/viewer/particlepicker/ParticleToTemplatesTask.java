package xmipp.viewer.particlepicker;

import xmipp.utils.Task;
import xmipp.viewer.particlepicker.training.gui.TemplatesJDialog;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;

public class ParticleToTemplatesTask implements Task
{

	private static TemplatesJDialog dialog;
	private TrainingParticle particle;

	
	public ParticleToTemplatesTask(TrainingParticle particle)
	{
		this.particle = particle;
	}
	
	public static void setTemplatesDialog(TemplatesJDialog d)
	{
		dialog = d;
	}
	
	@Override
	public void doTask()
	{
		try
		{
			Family family = particle.getFamily();
			family.addParticleToTemplates(particle);

			
			if(dialog != null && dialog.isVisible())
				dialog.loadTemplates(true);
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}

	}

}
