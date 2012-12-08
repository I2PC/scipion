package xmipp.particlepicker.tiltpair.gui;

import ij.gui.ImageWindow;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseWheelListener;
import java.util.List;
import javax.swing.SwingUtilities;

import xmipp.particlepicker.Micrograph;
import xmipp.particlepicker.ParticlePickerCanvas;
import xmipp.particlepicker.ParticlePickerJFrame;
import xmipp.particlepicker.tiltpair.model.TiltedMicrograph;
import xmipp.particlepicker.tiltpair.model.TiltedParticle;
import xmipp.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.particlepicker.tiltpair.model.UntiltedParticle;
import xmipp.particlepicker.training.model.TrainingMicrograph;
import xmipp.particlepicker.training.model.TrainingParticle;
import xmipp.utils.XmippWindowUtil;
import xmipp.jni.Particle;

public class TiltedMicrographCanvas extends ParticlePickerCanvas
{

	private TiltPairPickerJFrame frame;
	private UntiltedMicrograph um;
	private UntiltedMicrographCanvas uc;
	private TiltedParticle active;
	private boolean reload;
	private boolean drawalpha;

	public TiltedMicrographCanvas(TiltPairPickerJFrame frame)
	{
		super(frame.getMicrograph().getTiltedMicrograph().getImagePlus(frame.getParticlePicker().getFilters()));
		this.um = frame.getMicrograph();
		this.frame = frame;
		this.uc = (UntiltedMicrographCanvas)frame.getCanvas();
		//XmippWindowUtil.setLocation(0.7f, 0, iw);
		um.getTiltedMicrograph().runImageJFilters(frame.getParticlePicker().getFilters());
	}

	public void updateMicrograph()
	{
		this.um = frame.getMicrograph();
		updateMicrographData();
		active = null;
	}

	/**
	 * Adds particle or updates its position if onpick. If ondeletepick removes
	 * particle. Considers owner for selection to the first particle containing
	 * point. Sets dragged if onpick
	 */

	public void mousePressed(int x, int y)
	{
		setupScroll(x, y);
	}

	/**
	 * Updates particle position and repaints. Sets dragged to null at the end
	 */
	public void mouseReleased(MouseEvent e)
	{
		super.mouseReleased(e);
		if (reload)
			um.initAligner();
		reload = false;
	}

	/**
	 * Updates particle position and repaints if onpick.
	 */
	public void mouseDragged(int x, int y)
	{
		scroll(x, y);
	}

	public void mouseWheelMoved(int x, int y, int rotation)
	{

		if (rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);
		if (getMagnification() <= 1.0)
			imp.repaintWindow();

	}

	public void paint(Graphics g)
	{
		Graphics offgc;
		Image offscreen = null;
		Dimension d = getSize();

		// create the offscreen buffer and associated Graphics
		offscreen = createImage(d.width, d.height);
		offgc = offscreen.getGraphics();
		super.paint(offgc);
		Graphics2D g2 = (Graphics2D) offgc;
		g2.setColor(frame.getColor());
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
		if(frame.drawAngles())
			drawLine(Math.toRadians(um.getTiltedAngle()), g2);
		g.drawImage(offscreen, 0, 0, this);
	}

	public void mousePressed(MouseEvent e)
	{
		super.mousePressed(e);
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (frame.isPickingAvailable(e))
		{
			TiltedParticle p = um.getTiltedMicrograph().getParticle(x, y, (int) (frame.getParticleSize()));
			if (p != null)
			{
				if (SwingUtilities.isLeftMouseButton(e) && e.isShiftDown())
				{
					um.removeParticle(p.getUntiltedParticle());
					frame.updateMicrographsModel();
					if (p.getUntiltedParticle().isAdded())
						reload = true;
					frame.getCanvas().repaint();
				}
				else if (SwingUtilities.isLeftMouseButton(e))
					active = p;
			}
			else if (uc.hasActiveParticle() && SwingUtilities.isLeftMouseButton(e) && um.fits(x, y, frame.getParticleSize()))
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
	
	

	@Override
	public void mouseDragged(MouseEvent e)
	{

		super.mouseDragged(e);
		if (SwingUtilities.isRightMouseButton(e))
			return;//dont move particle if it is dragging image
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (active != null && um.fits(x, y, frame.getParticleSize()))
		{
			moveActiveParticle(x, y);
			if (active.getUntiltedParticle().isAdded())
				reload = true;
		}
		frame.setChanged(true);
		repaint();
	}



	@Override
	public void setActive(TrainingParticle p)
	{
		frame.getCanvas().setActive(((TiltedParticle) p).getUntiltedParticle());
	}

	@Override
	public ParticlePickerJFrame getFrame()
	{
		return frame;
	}

	@Override
	public Micrograph getMicrograph()
	{
		return um.getTiltedMicrograph();
	}

	@Override
	public TrainingParticle getActive()
	{
		return active;
	}
	


}
