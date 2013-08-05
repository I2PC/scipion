package xmipp.viewer.particlepicker;

import xmipp.utils.Task;
import xmipp.viewer.particlepicker.training.gui.TemplatesJDialog;
import xmipp.viewer.particlepicker.training.model.SingleParticlePicker;


public class UpdateTemplatesTask implements Task
{

	private static TemplatesJDialog dialog;
	private SingleParticlePicker picker;
	private int num;

	public UpdateTemplatesTask(SingleParticlePicker picker, int num)
	{
		this.picker = picker;
		this.num = num;
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
			picker.updateTemplates(num);
//			if(dialog != null && dialog.isVisible())
//				dialog.loadTemplates(true);
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}

	}

}
