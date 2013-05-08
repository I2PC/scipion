package xmipp.viewer.particlepicker;


import xmipp.jni.ImageGeneric;
import xmipp.utils.Task;
import xmipp.utils.TasksManager;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;

public class ParticleOutOfTemplatesTask implements Task
{
	
	private TrainingParticle particle;

	public ParticleOutOfTemplatesTask(TrainingParticle p)
	{
		this.particle = p;
	}

	@Override
	public void doTask()
	{
		
		try
		{
			
			Family family = particle.getFamily();
			
			ImageGeneric igp = particle.getImageGeneric();
//			System.out.printf("removing: %d %.2f %.2f %.2f\n", particle.getTemplateIndex(), particle.getTemplateRotation(), particle.getTemplateTilt(), particle.getTemplatePsi());
			family.getTemplates().removeAlignment(igp, particle.getTemplateIndex(), particle.getTemplateRotation(), particle.getTemplateTilt(), particle.getTemplatePsi());
			family.saveTemplates();
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}

	}

}
