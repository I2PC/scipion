package xmipp.viewer.particlepicker;

import xmipp.utils.Task;
import xmipp.viewer.particlepicker.training.gui.TemplatesJDialog;
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

			f.updateTemplates(picker);
			if(dialog != null && dialog.isVisible())
				dialog.loadTemplates(true);
			

		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}

	}

}
