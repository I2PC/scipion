package xmipp.viewer.particlepicker;

import xmipp.jni.ImageGeneric;
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
			if(dialog != null && dialog.isVisible())
				dialog.loadTemplates(true);

		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}

	}

}
