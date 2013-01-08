package xmipp.particlepicker.extract;

import ij.ImagePlus;

import java.awt.Graphics2D;

import xmipp.jni.Particle;
import xmipp.particlepicker.Micrograph;
import xmipp.particlepicker.ParticlePickerCanvas;
import xmipp.particlepicker.ParticlePickerJFrame;
import xmipp.particlepicker.training.model.TrainingParticle;

public class ExtractCanvas extends ParticlePickerCanvas
{

	private ExtractPickerJFrame frame;
	private ExtractMicrograph micrograph;
	private ExtractParticle active;

	public ExtractCanvas(ExtractPickerJFrame frame)
	{
		super(frame.getMicrograph().getImagePlus());
		this.frame = frame;
		this.micrograph = frame.getMicrograph();
		// TODO Auto-generated constructor stub
	}

	@Override
	public void refreshActive(Particle p)
	{
		active = (ExtractParticle)p;
		repaint();
		
	}

	@Override
	public TrainingParticle getActive()
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ParticlePickerJFrame getFrame()
	{
		return frame;
	}

	@Override
	public Micrograph getMicrograph()
	{
		return micrograph;
	}

	@Override
	protected void doCustomPaint(Graphics2D g2)
	{
		g2.setColor(frame.getColor());

		for (ExtractParticle p : micrograph.getParticles())
			drawShape(g2, p.getX(), p.getY(), frame.getParticlePicker().getSize(), false);
		
	}

	@Override
	public void updateMicrograph()
	{
		this.micrograph = frame.getMicrograph();
		updateMicrographData();
	}

}
