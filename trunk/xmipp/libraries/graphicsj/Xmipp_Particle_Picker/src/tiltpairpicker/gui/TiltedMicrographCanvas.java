package tiltpairpicker.gui;

import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.List;

import javax.swing.SwingUtilities;

import tiltpairpicker.model.TiltPairPicker;
import tiltpairpicker.model.TiltedMicrograph;
import tiltpairpicker.model.TiltedParticle;
import tiltpairpicker.model.UntiltedMicrograph;
import tiltpairpicker.model.UntiltedParticle;
import trainingpicker.gui.WindowUtils;
import trainingpicker.model.MicrographParticle;
import trainingpicker.model.TrainingParticle;

public class TiltedMicrographCanvas extends ImageCanvas implements MouseListener
{

	private TiltPairPickerJFrame frame;
	private UntiltedMicrograph untiltedmicrograph;
	private TiltedParticle dragged;
	private ImageWindow iw;
	private boolean reload;

	public TiltedMicrographCanvas(TiltPairPickerJFrame frame)
	{
		super(frame.getUntiltedMicrograph().getTiltedMicrograph().getImage());
		this.untiltedmicrograph = frame.getUntiltedMicrograph();
		this.frame = frame;
		iw = new ImageWindow(imp, this);
		WindowUtils.centerScreen(0.7, iw);

	}

	public void updateMicrograph()
	{
		this.untiltedmicrograph = frame.getUntiltedMicrograph();
		iw.setImage(untiltedmicrograph.getTiltedMicrograph().getImage());
		iw.updateImage(untiltedmicrograph.getTiltedMicrograph().getImage());
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
			untiltedmicrograph.initAligner();
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
		List<TiltedParticle> particles = untiltedmicrograph.getTiltedMicrograph().getParticles();
		for (TiltedParticle p : particles)
		{
			drawShape(g2, p, x0, y0, index == (particles.size() - 1));
			index++;
		}
		if (untiltedmicrograph.getActiveTiltedParticle() != null)
		{
			g2.setColor(Color.red);
			drawShape(g2, untiltedmicrograph.getActiveParticle().getTiltedParticle(), x0, y0, true);
		}
	}

	private void drawShape(Graphics2D g2, MicrographParticle p, int x0, int y0, boolean all)
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
		TiltedParticle p = untiltedmicrograph.getTiltedMicrograph().getParticle(x, y, (int) (frame.getParticleSize()));
		if (p != null)
		{
			if (SwingUtilities.isLeftMouseButton(e) && e.isControlDown())
			{
				untiltedmicrograph.getTiltedMicrograph().removeParticle(p);
				frame.updateMicrographsModel();
				if(p.getUntiltedParticle().isAdded())
					reload = true;
			}
			else if (SwingUtilities.isLeftMouseButton(e))
				dragged = p;
		}
		else if (untiltedmicrograph.hasActiveParticle() && SwingUtilities.isLeftMouseButton(e) && MicrographParticle.boxContainedOnImage(x, y, frame.getParticleSize(), imp))
		{
			UntiltedParticle activeparticle = untiltedmicrograph.getActiveParticle();
			if (activeparticle.getTiltedParticle() != null)
			{
				p = activeparticle.getTiltedParticle();
				p.setX(x);
				p.setY(y);
			}
			else
			{
				p = new TiltedParticle(x, y, untiltedmicrograph.getActiveParticle());

				untiltedmicrograph.getActiveParticle().setTiltedParticle(p);
				untiltedmicrograph.getTiltedMicrograph().addParticle(p);
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
		if (dragged != null && MicrographParticle.boxContainedOnImage(x, y, frame.getParticleSize(), imp))
		{
			dragged.setPosition(x, y);
			if(dragged.getUntiltedParticle().isAdded())
				reload = true;
		}
		frame.setChanged(true);
		repaint();
	}

}
