package xmipp.viewer.particlepicker.training.model;

import javax.swing.SwingWorker;

import xmipp.viewer.particlepicker.training.gui.SupervisedParticlePickerCanvas;

public class CenterParticleTask extends SwingWorker<String, Object>
{

	private static SupervisedParticlePickerCanvas canvas;
	private ManualParticle particle;
        private SupervisedParticlePicker ppicker;

	
	public CenterParticleTask(SupervisedParticlePickerCanvas canvas, SupervisedParticlePicker ppicker, ManualParticle particle)
	{
            this.canvas = canvas;
            this.ppicker = ppicker;
            this.particle = particle;
	}
	
	
	
	@Override
	protected String doInBackground() throws Exception
	{
		try
		{
			
			ppicker.centerParticle(particle);
			
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}
		
		return "";
	}
	
	 @Override
     protected void done() {
		 canvas.repaint();
     }

}
