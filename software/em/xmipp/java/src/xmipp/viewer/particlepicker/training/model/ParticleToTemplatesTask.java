package xmipp.viewer.particlepicker.training.model;

import javax.swing.SwingWorker;
import xmipp.viewer.particlepicker.training.gui.TemplatesJDialog;

public class ParticleToTemplatesTask extends SwingWorker<String, Object>
{

	private static TemplatesJDialog dialog;
	private ManualParticle particle;

	
	public ParticleToTemplatesTask(ManualParticle particle)
	{
		this.particle = particle;
	}
	
	public static void setTemplatesDialog(TemplatesJDialog d)
	{
		dialog = d;
	}
	
	@Override
	protected String doInBackground() throws Exception
	{
		try
		{
			
			SupervisedParticlePicker picker = (SupervisedParticlePicker)particle.getParticlePicker();
			if (particle.getMicrograph().fits(particle.getX(), particle.getY(), picker.getSize()))
				picker.addParticleToTemplates(particle);
			
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}
		
		return "";
	}
	
	 @Override
     protected void done() {
		 if(dialog != null && dialog.isVisible())
				dialog.loadTemplates(true);
     }

}
