package tiltpairpicker.gui;

import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

import javax.swing.SwingUtilities;

import tiltpairpicker.model.TiltPairPicker;
import tiltpairpicker.model.TiltedParticle;
import tiltpairpicker.model.UntiltedMicrograph;
import tiltpairpicker.model.UntiltedParticle;
import trainingpicker.gui.WindowUtils;
import trainingpicker.model.Constants;
import trainingpicker.model.MicrographParticle;
import xmipp.Particle;
import xmipp.TiltPairAligner;

public class UntiltedMicrographCanvas extends ImageCanvas implements MouseWheelListener
{

	private TiltPairPickerJFrame frame;
	private TiltedMicrographCanvas tiltedcanvas;
	private UntiltedParticle active;
	private TiltPairPicker pppicker;
	private UntiltedMicrograph untiltedmicrograph;
	private ImageWindow iw;
	private boolean reload = false;
	

	public UntiltedMicrographCanvas(TiltPairPickerJFrame frame)
	{
		super(frame.getMicrograph().getImage());
		this.untiltedmicrograph = frame.getMicrograph();
	
		this.frame = frame;
		addMouseWheelListener(this);
		this.pppicker = frame.getParticlePairPicker();
		tiltedcanvas = new TiltedMicrographCanvas(frame);
		iw = new ImageWindow(imp, this);
		WindowUtils.centerScreen(0, iw);
	}

	public void updateMicrograph()
	{
		this.untiltedmicrograph = frame.getMicrograph();
		iw.setImage(untiltedmicrograph.getImage());
		iw.updateImage(untiltedmicrograph.getImage());
		tiltedcanvas.updateMicrograph();
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
		{
			setupScroll(x, y);
			tiltedcanvas.mousePressed(x, y);
		}
		if(active != null && !active.isAdded() && active.getTiltedParticle() != null)
				untiltedmicrograph.addParticleToAligner(active);
		UntiltedParticle p = untiltedmicrograph.getParticle(x, y, (int) (frame.getParticleSize()));
		
		if (p != null)
		{
			if (SwingUtilities.isLeftMouseButton(e) && e.isControlDown())
			{
				untiltedmicrograph.removeParticle(p);
				frame.updateMicrographsModel();
				if(p.isAdded())
					untiltedmicrograph.initAligner();
			}
			else if (SwingUtilities.isLeftMouseButton(e))
				setActive(p);
		}
		else if (SwingUtilities.isLeftMouseButton(e) && MicrographParticle.boxContainedOnImage(x, y, frame.getParticleSize(), imp))
		{
			p = new UntiltedParticle(x, y, untiltedmicrograph);
			untiltedmicrograph.addParticle(p);
			setActive(p);
			frame.updateMicrographsModel();
		}
		frame.setChanged(true);
		repaint();
		tiltedcanvas.repaint();

	}
	
	public void setActive(UntiltedParticle up)
	{
		active = up;
		untiltedmicrograph.setActiveParticle(active);
		if(active.getTiltedParticle() == null && untiltedmicrograph.getAddedCount() >= 4)
			untiltedmicrograph.setAlignerTiltedParticle(up);
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
			tiltedcanvas.mouseDragged(e.getX(), e.getY());
			return;
		}
		if (active != null && MicrographParticle.boxContainedOnImage(x, y, frame.getParticleSize(), imp))
		{
			active.setPosition(x, y);
			if(active.isAdded())
				reload = true;
		}
		frame.setChanged(true);
		repaint();

	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent e)
	{
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		tiltedcanvas.setMagnification(magnification);
		int rotation = e.getWheelRotation();
		if (rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);
		
		tiltedcanvas.mouseWheelMoved(x, y, rotation);
		
	}

	public TiltedMicrographCanvas getTiltedCanvas()
	{
		return tiltedcanvas;
	}

	public void paint(Graphics g)
	{
		super.paint(g);
		Graphics2D g2 = (Graphics2D) g;
		g2.setColor(frame.getColor());
		int x0 = (int) getSrcRect().getX();
		int y0 = (int) getSrcRect().getY();
		int index = 0;

		for (MicrographParticle p : untiltedmicrograph.getParticles())
		{
			drawShape(g2, p, x0, y0, index == (untiltedmicrograph.getParticles().size() - 1));
			index++;
		}
		if (untiltedmicrograph.getActiveParticle() != null)
		{
			g2.setColor(Color.red);
			drawShape(g2, untiltedmicrograph.getActiveParticle(), x0, y0, true);
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
	


}
