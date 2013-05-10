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
			Family family = particle.getFamily();

			ImageGeneric igp = particle.getImageGeneric();
			// will happen only in manual mode
			if (family.getTemplateIndex() < family.getTemplatesNumber())
				family.setTemplate(igp);
			else
			{
				
				double[] align = family.getTemplates().alignImage(igp);
				family.applyAlignment(particle, igp, align);

			}
			family.saveTemplates();
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}

	}

}
