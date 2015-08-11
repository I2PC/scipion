package xmipp.viewer.particlepicker.tiltpair.gui;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.util.List;

import javax.swing.SwingUtilities;

import xmipp.jni.Particle;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.tiltpair.model.TiltPairPicker;
import xmipp.viewer.particlepicker.tiltpair.model.TiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.TiltedParticle;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedParticle;
import xmipp.viewer.particlepicker.training.model.ManualParticle;

public class TiltedMicrographCanvas extends ParticlePickerCanvas
{

	private UntiltedMicrograph um;
	private UntiltedMicrographCanvas uc;
	private TiltedParticle active;

	public TiltedMicrographCanvas(TiltPairPickerJFrame frame)
	{
		super(frame, frame.getMicrograph().getTiltedMicrograph());
		this.um = frame.getMicrograph();
		this.uc = (UntiltedMicrographCanvas) frame.getCanvas();
		//XmippWindowUtil.setLocation(0.7f, 0, iw);
	}
	
	

	public void updateMicrograph()
	{
		um = getFrame().getMicrograph();
		TiltedMicrograph m = getFrame().getMicrograph().getTiltedMicrograph();
		updateMicrograph(m);

	}

	/**
	 * Adds particle or updates its position if onpick. If ondeletepick removes
	 * particle. Considers owner for selection to the first particle containing
	 * point. Sets dragged if onpick
	 */

	

	public void mousePressed(MouseEvent e)
	{
		super.mousePressed(e);
		if (frame.isPickingAvailable(e))
		{
			int x = super.offScreenX(e.getX());
			int y = super.offScreenY(e.getY());

			if (frame.isPickingAvailable(e))
			{
				TiltedParticle p = um.getTiltedMicrograph().getParticle(x, y, (int) (getFrame().getParticleSize()));
				if (p != null)
				{
					if (SwingUtilities.isLeftMouseButton(e) && e.isShiftDown())
					{
						um.removeParticle(p.getUntiltedParticle());
						frame.updateMicrographsModel();
						frame.getCanvas().repaint();
					}
					else if (SwingUtilities.isLeftMouseButton(e))
						active = p;
				}
				else if (uc.hasActiveParticle() && SwingUtilities.isLeftMouseButton(e) && um.fits(x, y, getFrame().getParticleSize()))
				{
					UntiltedParticle uactive = uc.getActiveParticle();
					if (uactive.getTiltedParticle() != null)
						p = uactive.getTiltedParticle();
					else
					{
						p = new TiltedParticle(x, y, uc.getActiveParticle());

						uc.getActiveParticle().setTiltedParticle(p);
						um.getTiltedMicrograph().addParticle(p);
					}
					active = p;
					frame.updateMicrographsModel();
				}
				frame.setChanged(true);
				repaint();
			}
		}
	}

	public void mousePressed(int x, int y)
	{
		setupScroll(x, y);
	}

	/**
	 * Updates particle position and repaints if onpick.
	 */
	public void mouseDragged(int x, int y)
	{
		scroll(x, y);
	}

	@Override
	public void mouseDragged(MouseEvent e)
	{

		super.mouseDragged(e);
		if (frame.isPickingAvailable(e))
		{
			int x = super.offScreenX(e.getX());
			int y = super.offScreenY(e.getY());

			if (active != null && um.fits(x, y, getFrame().getParticleSize()))
			{
				setActiveMoved(true);
				moveActiveParticle(x, y);

			}
			frame.setChanged(true);
			repaint();
		}
	}

	public void mouseReleased(MouseEvent e)
	{
		super.mouseReleased(e);
		if (frame.isPickingAvailable(e))
		{
			super.mouseReleased(e);
			int x = super.offScreenX(e.getX());
			int y = super.offScreenY(e.getY());
			manageActive(x, y);
		}
	}
        
        @Override
	public void mouseWheelMoved(MouseWheelEvent e)
	{
		super.mouseWheelMoved(e);
		if (!e.isShiftDown())
			return;
		if(getParticlePicker().getZoom() != uc.getMagnification())
			uc.mouseWheelMoved(e);	
	}

	@Override
	public void refreshActive(Particle p)
	{
		if (p != null)
			frame.getCanvas().refreshActive(((TiltedParticle) p).getUntiltedParticle());
	}

	@Override
	public TiltPairPickerJFrame getFrame()
	{
		return (TiltPairPickerJFrame)frame;
	}

	@Override
	public TiltedMicrograph getMicrograph()
	{
		return (TiltedMicrograph)micrograph;
	}

	@Override
	public ManualParticle getActive()
	{
		return active;
	}

	@Override
	protected void doCustomPaint(Graphics2D g2)
	{
		g2.setColor(getFrame().getColor());
		int index = 0;
		List<TiltedParticle> particles = um.getTiltedMicrograph().getParticles();
		for (TiltedParticle p : particles)
		{
			drawShape(g2, p, index == (particles.size() - 1));
			index++;
		}

		if (uc.getActiveTiltedParticle() != null)
		{
			g2.setColor(Color.red);
			drawShape(g2, uc.getActiveParticle().getTiltedParticle(), true);
		}
		if (getFrame().drawAngles())
			drawLine(Math.toRadians(um.getTiltedAngle()), g2);

	}

	@Override
	protected Particle getLastParticle()
	{
		return uc.getLastParticle().getTiltedParticle();
	}

	protected void manageActive(int x, int y)
	{
		if (!activemoved)
			return;
		if (active.getUntiltedParticle().isAdded())
			um.initAligner();
		setActiveMoved(false);
	}

	

	@Override
	public TiltPairPicker getParticlePicker()
	{
		// TODO Auto-generated method stub
		return (TiltPairPicker)picker;
	}

}
