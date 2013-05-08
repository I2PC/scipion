package xmipp.viewer.particlepicker;

import xmipp.jni.ImageGeneric;
import xmipp.jni.Particle;
import xmipp.utils.TasksManager;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;

public class ParticleToTemplates implements Runnable
{

	private TrainingParticle particle;
	private int templateIndex;
	private boolean center;
	
	public ParticleToTemplates(TrainingParticle particle, boolean center)
	{
		this.particle = particle;
		this.center = center;
	}
	
	@Override
	public void run()
	{
		try
		{
			Particle shift = null;
			Family family = particle.getFamily();

			ImageGeneric igp = particle.getImageGeneric();
			// will happen only in manual mode
			if (family.getTemplateIndex() < family.getTemplatesNumber())
				family.setTemplate((int) (ImageGeneric.FIRST_IMAGE + templateIndex), igp);
			else
			{
				if (center)
				{
					shift = family.getTemplates().bestShift(igp);
					particle.setX(particle.getX() + shift.getX());
					particle.setY(particle.getY() + shift.getY());
				}
				double[] align = family.getTemplates().alignImage(igp);
				particle.setLastalign(align);
				System.out.printf("adding particle: %d %.2f %.2f %.2f\n", particle.getTemplateIndex(), particle.getTemplateRotation(), particle.getTemplateTilt(), particle.getTemplatePsi());

			}
			family.saveTemplates();
			TasksManager.getInstance().removeTask(this);
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}

	}

}
