package xmipp.viewer.particlepicker;

import xmipp.utils.Task;
import xmipp.viewer.particlepicker.training.gui.TemplatesJDialog;


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
		try
		{
			picker.updateTemplates();
			if(dialog != null && dialog.isVisible())
				dialog.loadTemplates(true);
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}

	}

}
