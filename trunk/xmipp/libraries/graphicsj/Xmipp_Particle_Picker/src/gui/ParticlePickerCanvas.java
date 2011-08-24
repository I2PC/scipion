package gui;

import ij.gui.ImageCanvas;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

import javax.swing.SwingUtilities;

import model.Micrograph;
import model.MicrographFamilyData;
import model.Particle;

public class ParticlePickerCanvas extends ImageCanvas implements
		MouseWheelListener {

	private ParticlePickerJFrame frame;
	private Micrograph micrograph;
	private Particle dragged;
	final static BasicStroke dashedst = new BasicStroke(1.0f,
			BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f,
			new float[] { 10.0f }, 0.0f);
	final static BasicStroke continuousst = new BasicStroke();

	public ParticlePickerCanvas(ParticlePickerJFrame frame) {
		super(frame.getMicrograph().getImage());
		this.micrograph = frame.getMicrograph();
		this.frame = frame;
		addMouseWheelListener(this);

		// TODO Auto-generated constructor stub
	}

	/**
	 * Adds particle or updates its position if onpick. If ondeletepick removes
	 * particle. Considers owner for selection to the first particle containing
	 * point. Sets dragged if onpick
	 */
	public void mousePressed(MouseEvent e) {
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if (frame.getTool() != Tool.PICKER) {
			super.mousePressed(e);
			return;
		}
		if (SwingUtilities.isRightMouseButton(e)) {
			setupScroll(x, y);
			return;
		}
		Particle p = micrograph.getParticle(x, y);
		if (p != null) {
			if (SwingUtilities.isLeftMouseButton(e) && e.isControlDown()) {
				micrograph.removeParticle(p);
				frame.updateMicrographsModel();
			} else if (SwingUtilities.isLeftMouseButton(e)) {
				p.setPosition(x, y);
				dragged = p;
			}
		} else if (SwingUtilities.isLeftMouseButton(e)
				&& Particle.boxContainedOnImage(x, y, frame.getFamily()
						.getSize(), imp)) {
			p = new Particle(x, y, frame.getFamily(), micrograph);
			micrograph.addParticle(p);
			dragged = p;
			frame.updateMicrographsModel();
		}
		frame.setChanged(true);
		repaint();
	}

	/**
	 * Updates particle position and repaints. Sets dragged to null at the end
	 */
	public void mouseReleased(MouseEvent e) {
		if (frame.getTool() != Tool.PICKER) {
			super.mouseReleased(e);
			return;
		}

		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if (dragged == null)// not onpick
			return;
		Particle p = null;
		p = dragged;
		dragged = null;
		if (!Particle.boxContainedOnImage(x, y, p.getFamily().getSize(), imp))
			return;
		p.setPosition(x, y);
		frame.setChanged(true);
		repaint();
	}

	/**
	 * Updates particle position and repaints if onpick.
	 */
	@Override
	public void mouseDragged(MouseEvent e) {
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if (frame.getTool() != Tool.PICKER) {
			super.mouseDragged(e);
			return;
		}
		if (SwingUtilities.isRightMouseButton(e)) {
			scroll(e.getX(), e.getY());
			return;
		}
		if (dragged == null)
			return;

		if (!Particle.boxContainedOnImage(x, y, dragged.getFamily().getSize(),
				imp))
			return;
		dragged.setPosition(x, y);
		frame.setChanged(true);
		repaint();
	}

	public void paint(Graphics g) {
		super.paint(g);
		Graphics2D g2 = (Graphics2D) g;
		int x0 = (int) getSrcRect().getX();
		int y0 = (int) getSrcRect().getY();
		int radius;
		int count = 0;
		int x, y;

		for (MicrographFamilyData mfdata : micrograph.getFamiliesData()) {
			g2.setColor(mfdata.getFamily().getColor());
			radius = (int) (mfdata.getFamily().getSize() / 2 * magnification);
			for (Particle p : mfdata.getParticles()) {

				if (p.isAuto())
					g2.setStroke(dashedst);
				else
					g2.setStroke(continuousst);

				count++;
				x = (int) ((p.getX() - x0) * magnification);
				y = (int) ((p.getY() - y0) * magnification);
				
				drawShape(g2, x, y, radius, count);
			}
		}
	}
	
	void drawShape(Graphics2D g2, int x, int y, int radius, int label)
	{
		if (frame.isShapeSelected(Shape.Rectangle))
			g2.drawRect(x - radius, y - radius, radius * 2, radius * 2);
		if (frame.isShapeSelected(Shape.Circle))
			g2.drawOval(x - radius, y - radius, radius * 2, radius * 2);
		if (frame.isShapeSelected(Shape.Center))
		{
			g2.drawLine(x - 2, y - 2, x + 2, y + 2);
			g2.drawLine(x + 2, y - 2, x - 2, y + 2);
		}
		g2.drawString(Integer.toString(label), x, y - radius);
	}

	public void updateMicrograph() {
		this.micrograph = frame.getMicrograph();
		imp = micrograph.getImage();
		setImageUpdated();
		repaint();
	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent e) {
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		int rotation = e.getWheelRotation();
		if (rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);

	}

}
