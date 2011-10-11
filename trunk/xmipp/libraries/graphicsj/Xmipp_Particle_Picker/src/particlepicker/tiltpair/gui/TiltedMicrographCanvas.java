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
	private TiltedParticle dragged;
	private ImageWindow iw;
	private boolean reload;
	

	public TiltedMicrographCanvas(TiltPairPickerJFrame frame)
	{
		super(frame.getMicrograph().getTiltedMicrograph().getImage());
		this.um = frame.getMicrograph();
		this.frame = frame;
		iw = new ImageWindow(imp, this);
		WindowUtils.centerScreen(0.7, 0, iw);
		addMouseWheelListener(this);

	}

	public void updateMicrograph()
	{
		this.um = frame.getMicrograph();
		iw.setImage(um.getTiltedMicrograph().getImage());
		iw.updateImage(um.getTiltedMicrograph().getImage());
		dragged = null;
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
		if (frame.getTool() != Tool.PICKER)
		{
			super.mouseReleased(e);
			return;
		}
		if(reload)
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
			drawShape(g2, p, x0, y0, index == (particles.size() - 1));
			index++;
		}
		if (um.getActiveTiltedParticle() != null)
		{
			g2.setColor(Color.red);
			drawShape(g2, um.getActiveParticle().getTiltedParticle(), x0, y0, true);
		}
	}

	private void drawShape(Graphics2D g2, TrainingParticle p, int x0, int y0, boolean all)
	{
		int size = (int) (frame.getParticleSize() * magnification);
		int radius = (int) (frame.getParticleSize() / 2 * magnification);
		int x = (int) ((p.getX() - x0) * magnification);
		int y = (int) ((p.getY() - y0) * magnification);
		int distance = (int) (5 * magnification);

		if (frame.isShapeSelected(Shape.Rectangle) || all)
			g2.drawRect(x - radius, y - radius, size, size);
		if (frame.isShapeSelected(Shape.Circle) || all)
			g2.drawOval(x - radius, y - radius, size, size);
		if (frame.isShapeSelected(Shape.Center) || all)
		{
			g2.drawLine(x, y - distance, x, y + distance);
			g2.drawLine(x + distance, y, x - distance, y);
		}
	}

	public void mousePressed(MouseEvent e)
	{
		if (frame.getTool() != Tool.PICKER)
		{
			super.mousePressed(e);
			return;
		}
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (SwingUtilities.isRightMouseButton(e))
			setupScroll(x, y);
		TiltedParticle p = um.getTiltedMicrograph().getParticle(x, y, (int) (frame.getParticleSize()));
		if (p != null)
		{
			if (SwingUtilities.isLeftMouseButton(e) && e.isControlDown())
			{
				um.removeParticle(p.getUntiltedParticle());
				frame.updateMicrographsModel();
				if(p.getUntiltedParticle().isAdded())
					reload = true;
				frame.getCanvas().repaint();
			}
			else if (SwingUtilities.isLeftMouseButton(e))
				dragged = p;
		}
		else if (um.hasActiveParticle() && SwingUtilities.isLeftMouseButton(e) && Particle.boxContainedOnImage(x, y, frame.getParticleSize(), imp))
		{
			UntiltedParticle active = um.getActiveParticle();
			if (active.getTiltedParticle() != null)
			{
				p = active.getTiltedParticle();
//				p.setX(x);
//				p.setY(y);
			}
			else
			{
				p = new TiltedParticle(x, y, um.getActiveParticle());

				um.getActiveParticle().setTiltedParticle(p);
				um.getTiltedMicrograph().addParticle(p);
			}
			dragged = p;
			frame.updateMicrographsModel();
		}
		frame.setChanged(true);
		repaint();

	}

	@Override
	public void mouseDragged(MouseEvent e)
	{

		if (frame.getTool() != Tool.PICKER)
		{
			super.mouseDragged(e);
			return;
		}
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if (SwingUtilities.isRightMouseButton(e))
		{
			scroll(e.getX(), e.getY());
			return;
		}
		if (dragged != null && Particle.boxContainedOnImage(x, y, frame.getParticleSize(), imp))
		{
			dragged.setPosition(x, y);
			if(dragged.getUntiltedParticle().isAdded())
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
		if(x0 < 0)
			x0 = 0;
		if(x0 + width > imp.getWidth())
			x0 = imp.getWidth() - width;
		int y0 = p.getY() - height / 2;
		if(y0 < 0)
			y0 = 0;
		if(y0 + height > imp.getHeight())
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
	public void mouseWheelMoved(MouseWheelEvent e)
	{
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		int rotation = e.getWheelRotation();
		if (rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);
		
	}

}
