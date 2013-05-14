package xmipp.viewer.particlepicker;

import java.util.List;

import xmipp.jni.ImageGeneric;
import xmipp.utils.Task;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.TrainingMicrograph;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;

public class UpdateTemplatesTask implements Task
{

	private SingleParticlePicker picker;

	public UpdateTemplatesTask(SingleParticlePicker picker)
	{
		this.picker = picker;
	}

	@Override
	public void doTask()
	{
		if (picker.getMode() != Mode.Manual)
			return;

		if (!picker.hasManualParticles())
			return;

		picker.initTemplates();
		ImageGeneric igp;
		List<TrainingParticle> particles;
		TrainingParticle particle;
		double[] align;
		try
		{
			for (TrainingMicrograph m : picker.getMicrographs())
			{
				for (int i = 0; i < m.getManualParticles().size(); i++)
				{
					particles = m.getManualParticles();
					particle = particles.get(i);
					igp = particle.getImageGeneric();
					if (picker.getTemplateIndex() < picker.getTemplatesNumber())
						picker.setTemplate(igp);
					else
					{
						align = picker.getTemplates().alignImage(igp);
						picker.applyAlignment(particle, igp, align);
					}
				}
			}
			
			picker.saveTemplates();
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}


	}

}
