package xmipp.viewer.particlepicker;

import xmipp.utils.Task;
import xmipp.viewer.particlepicker.training.gui.TemplatesJDialog;
import xmipp.viewer.particlepicker.training.model.ManualParticle;

public class ParticleToTemplatesTask implements Task
{

	private static TemplatesJDialog dialog;
	private ManualParticle particle;

	
	public ParticleToTemplatesTask(ManualParticle particle)
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
			
			SingleParticlePicker picker = (SingleParticlePicker)particle.getParticlePicker();

			picker.addParticleToTemplates(particle);
			if(dialog != null && dialog.isVisible())
				dialog.loadTemplates(true);
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}

	}

}
