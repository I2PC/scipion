package xmipp.viewer.particlepicker.training.model;

import javax.swing.SwingWorker;

import xmipp.viewer.particlepicker.training.gui.SupervisedPickerCanvas;

public class CenterParticleTask extends SwingWorker<String, Object>
{

	private static SupervisedPickerCanvas canvas;
	private ManualParticle particle;
        private SupervisedParticlePicker ppicker;

	
	public CenterParticleTask(SupervisedPickerCanvas canvas, SupervisedParticlePicker ppicker, ManualParticle particle)
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
