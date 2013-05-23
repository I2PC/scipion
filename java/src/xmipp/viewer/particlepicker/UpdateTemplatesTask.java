package xmipp.viewer.particlepicker;

import java.util.List;

import xmipp.jni.ImageGeneric;
import xmipp.jni.Particle;
import xmipp.utils.Task;
import xmipp.utils.TasksManager;
import xmipp.viewer.particlepicker.training.gui.TemplatesJDialog;
import xmipp.viewer.particlepicker.training.model.FamilyState;
import xmipp.viewer.particlepicker.training.model.MicrographFamilyData;
import xmipp.viewer.particlepicker.training.model.TrainingMicrograph;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;
import xmipp.viewer.particlepicker.training.model.TrainingPicker;

public class UpdateTemplatesTask implements Task
{

	private TrainingPicker picker;
	private static TemplatesJDialog dialog;

	public UpdateTemplatesTask(TrainingPicker picker)
	{
		this.picker = picker;
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
			Family f = picker.getFamily();

			if (f.getStep() != FamilyState.Manual)
				return;

			f.initTemplates();
			ImageGeneric igp;
			List<TrainingParticle> particles;
			MicrographFamilyData mfd;
			TrainingParticle particle;
			try
			{
				for (TrainingMicrograph m : picker.getMicrographs())
				{
					mfd = m.getFamilyData(f);
					for (int i = 0; i < mfd.getManualParticles().size(); i++)
					{
						particles = mfd.getManualParticles();
						particle = particles.get(i);
						igp = particle.getImageGeneric();
						if (f.getTemplateIndex() < f.getTemplatesNumber())
							f.setTemplate(igp);
						else
						{
							double[] align = f.getTemplates().alignImage(igp);
							particle.setLastalign(align);
						}
					}
				}
				f.saveTemplates();
				if(dialog != null && dialog.isVisible())
					dialog.loadTemplates(true);
			}
			catch (Exception e)
			{
				throw new IllegalArgumentException(e.getMessage());
			}

		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}

	}

}
