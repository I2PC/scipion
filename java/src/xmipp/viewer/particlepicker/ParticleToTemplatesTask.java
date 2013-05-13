package xmipp.viewer.particlepicker;

import xmipp.jni.ImageGeneric;
import xmipp.jni.Particle;
import xmipp.utils.Task;
import xmipp.utils.TasksManager;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;

public class ParticleToTemplatesTask implements Task
{

	private TrainingParticle particle;

	
	public ParticleToTemplatesTask(TrainingParticle particle)
	{
		this.particle = particle;
	}
	
	@Override
	public void doTask()
	{
		try
		{
			SingleParticlePicker picker = (SingleParticlePicker)particle.getParticlePicker();

			ImageGeneric igp = particle.getImageGeneric();
			// will happen only in manual mode
			if (picker.getTemplateIndex() < picker.getTemplatesNumber())
				picker.setTemplate(igp);
			else
			{
				
				double[] align = picker.getTemplates().alignImage(igp);
				picker.applyAlignment(particle, igp, align);

			}
			picker.saveTemplates();
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}

	}

}
