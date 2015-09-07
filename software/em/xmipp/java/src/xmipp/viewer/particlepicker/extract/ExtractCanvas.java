package xmipp.viewer.particlepicker.extract;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.event.MouseEvent;

import javax.swing.SwingUtilities;

import xmipp.jni.Particle;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.training.gui.SupervisedPickerJFrame;

public class ExtractCanvas extends ParticlePickerCanvas
{

	private ExtractParticle active;

	public ExtractCanvas(ExtractPickerJFrame frame)
	{
		super(frame);
		active = getLastParticle();
	}

	@Override
	public void refreshActive(Particle p)
	{
		active = (ExtractParticle) p;
		repaint();

	}

	@Override
	public ExtractParticle getActive()
	{
		return active;
	}



	@Override
	public ExtractMicrograph getMicrograph()
	{
		return (ExtractMicrograph)micrograph;
	}

	@Override
	protected void doCustomPaint(Graphics2D g2)
	{
		Color color;
		double score;
		for (ExtractParticle p : getMicrograph().getParticles())
		{
			score = p.getScore(getFrame().getColorHelper().getId());
			color = getFrame().getColorHelper().getColor(score);
			if (p.isEnabled())
				drawShape(g2, p.getX(), p.getY(), getFrame().getParticlePicker().getSize(), false, continuousst, color);
			else
				drawShape(g2, p.getX(), p.getY(), getFrame().getParticlePicker().getSize(), false, dashedst, color);
		}
		if (active != null)
		{
			score = active.getScore(getFrame().getColorHelper().getId());
			color = getFrame().getColorHelper().getColor(score);
			drawShape(g2, active.getX(), active.getY(), frame.getParticlePicker().getSize(), true, activest, color);
		}
	}

	@Override
	protected ExtractParticle getLastParticle()
	{
		if (getMicrograph().getParticles().isEmpty())
			return null;
		return getMicrograph().getParticles().get(getMicrograph().getParticles().size() - 1);
	}

	public void mousePressed(MouseEvent e)
	{
		super.mousePressed(e);
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (frame.isPickingAvailable(e))
		{
			if (frame.isEraserMode())
			{
				getMicrograph().removeParticles(x, y, picker);
				active = getLastParticle();
				refresh();

				return;
			}
			ExtractParticle p = getMicrograph().getParticle(x, y, picker.getSize());

			if (p != null)
			{
				if (SwingUtilities.isLeftMouseButton(e))
				{
					if (e.isShiftDown())
						p.setEnabled(!p.isEnabled());
					active = p;
				}
				getFrame().refreshActiveOnGallery(active);
			}

			refresh();
		}
	}

	@Override
	public void mouseDragged(MouseEvent e)
	{

		super.mouseDragged(e);
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if (frame.isPickingAvailable(e))
		{
			if (frame.isEraserMode())
			{
				getMicrograph().removeParticles(x, y, picker);
				active = getLastParticle();
				refresh();
				return;
			}
			if (active == null)
				return;

			if (!micrograph.fits(x, y, picker.getSize()))
				return;
			moveActiveParticle(x, y);
			active.setEnabled(true);// if it was disabled gets enabled
			getFrame().refreshActiveOnGallery(active);// if enabled propagate info
			frame.setChanged(true);
			repaint();
		}
	}

	
	
	protected void moveActiveParticle(int x, int y)
	{
		active.setEnabled(true);
		super.moveActiveParticle(x, y);
		getFrame().refreshActiveOnGallery(active);
	}

	@Override
	protected void manageActive(int x, int y)
	{
		if(!activemoved)
			return;
		if (!micrograph.fits(x, y, picker.getSize()))
			return;
		moveActiveParticle(x, y);
		setActiveMoved(false);
		
	}
	
	@Override
	public ExtractPickerJFrame getFrame()
	{
		return (ExtractPickerJFrame)frame;
	}

	@Override
	public ExtractParticlePicker getParticlePicker()
	{
		// TODO Auto-generated method stub
		return (ExtractParticlePicker)picker;
	}

}
