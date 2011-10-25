package particlepicker.tiltpair.gui;

import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.List;

import javax.swing.SwingUtilities;

import particlepicker.ParticlePickerCanvas;
import particlepicker.ParticlePickerJFrame;
import particlepicker.Tool;
import particlepicker.WindowUtils;
import particlepicker.tiltpair.model.TiltedParticle;
import particlepicker.tiltpair.model.UntiltedMicrograph;
import particlepicker.tiltpair.model.UntiltedParticle;
import particlepicker.training.model.TrainingParticle;
import xmipp.Particle;

public class TiltedMicrographCanvas extends ParticlePickerCanvas implements MouseListener, MouseWheelListener
{

	private TiltPairPickerJFrame frame;
	private UntiltedMicrograph um;
	private TiltedParticle active;
	private ImageWindow iw;
	private boolean reload;

	public TiltedMicrographCanvas(TiltPairPickerJFrame frame)
	{
		super(frame.getMicrograph().getTiltedMicrograph().getImagePlus());
		this.um = frame.getMicrograph();
		this.frame = frame;
		iw = new ImageWindow(imp, this);
		WindowUtils.centerScreen(0.7, 0, iw);
		addMouseWheelListener(this);

	}

	public void updateMicrograph()
	{
		this.um = frame.getMicrograph();
		iw.setImage(um.getTiltedMicrograph().getImagePlus());
		iw.updateImage(um.getTiltedMicrograph().getImagePlus());
		active = null;
	}

	public void mouseEntered(MouseEvent e)
	{
		if (frame.getTool() != Tool.PICKER)
		{
			super.mouseEntered(e);
			return;
		}
		setCursor(crosshairCursor);
	}

	public void mouseMoved(MouseEvent e)
	{
		if (frame.getTool() != Tool.PICKER)
		{
			super.mouseMoved(e);
			return;
		}
		setCursor(crosshairCursor);
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

	}

	public void paint(Graphics g)
	{
		super.paint(g);
		Graphics2D g2 = (Graphics2D) g;
		g2.setColor(frame.getColor());
		int x0 = (int) getSrcRect().getX();
		int y0 = (int) getSrcRect().getY();
		int index = 0;
		List<TiltedParticle> particles = um.getTiltedMicrograph().getParticles();
		for (TiltedParticle p : particles)
		{
			drawShape(g2, p, index == (particles.size() - 1));
			index++;
		}
		if (um.getActiveTiltedParticle() != null)
		{
			g2.setColor(Color.red);
			drawShape(g2, um.getActiveParticle().getTiltedParticle(), true);
		}
	}

	public void mousePressed(MouseEvent e)
	{
		super.mousePressed(e);
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (frame.isPickingAvailable())
		{
			TiltedParticle p = um.getTiltedMicrograph().getParticle(x, y, (int) (frame.getParticleSize()));
			if (p != null)
			{
				if (SwingUtilities.isLeftMouseButton(e) && e.isControlDown())
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
			else if (um.hasActiveParticle() && SwingUtilities.isLeftMouseButton(e) && Particle.boxContainedOnImage(x, y, frame.getParticleSize(), imp))
			{
				UntiltedParticle uactive = um.getActiveParticle();
				if (uactive.getTiltedParticle() != null)
					p = uactive.getTiltedParticle();
				else
				{
					p = new TiltedParticle(x, y, um.getActiveParticle());

					um.getActiveParticle().setTiltedParticle(p);
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
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (active != null && Particle.boxContainedOnImage(x, y, frame.getParticleSize(), imp))
		{
			active.setPosition(x, y);
			if (active.getUntiltedParticle().isAdded())
				reload = true;
		}
		frame.setChanged(true);
		repaint();
	}

	public void moveTo(TiltedParticle p)
	{

		int width = (int) getSrcRect().getWidth();
		int height = (int) getSrcRect().getHeight();
		int x0 = p.getX() - width / 2;
		if (x0 < 0)
			x0 = 0;
		if (x0 + width > imp.getWidth())
			x0 = imp.getWidth() - width;
		int y0 = p.getY() - height / 2;
		if (y0 < 0)
			y0 = 0;
		if (y0 + height > imp.getHeight())
			y0 = imp.getHeight() - height;
		Rectangle r = new Rectangle(x0, y0, width, height);
		setMagnification(frame.getCanvas().getMagnification());
		if (!getSrcRect().contains(r))
		{
			setSourceRect(r);
			repaint();
		}
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

}
