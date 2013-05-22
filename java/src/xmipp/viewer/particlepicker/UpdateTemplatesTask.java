package xmipp.viewer.particlepicker;

import java.util.List;

import xmipp.jni.ImageGeneric;
import xmipp.utils.Task;
import xmipp.viewer.particlepicker.training.gui.TemplatesJDialog;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.TrainingMicrograph;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;


public class UpdateTemplatesTask implements Task
{

	private static TemplatesJDialog dialog;
	private SingleParticlePicker picker;

	public UpdateTemplatesTask(SingleParticlePicker picker)
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
			picker.updateTemplates();
			if(dialog != null && dialog.isVisible())
				dialog.loadTemplates(true);
			System.out.println("Templates updated");
		



	}

}
