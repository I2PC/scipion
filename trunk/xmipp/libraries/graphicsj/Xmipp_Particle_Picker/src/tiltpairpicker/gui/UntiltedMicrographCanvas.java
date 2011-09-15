package tiltpairpicker.gui;

import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import java.awt.BasicStroke;
import java.awt.Cursor;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.List;
import javax.swing.SwingUtilities;

import tiltpairpicker.model.ParticlePairPicker;
import tiltpairpicker.model.UntiltedMicrograph;
import tiltpairpicker.model.UntiltedParticle;
import trainingpicker.gui.WindowUtils;
import trainingpicker.model.AutomaticParticle;
import trainingpicker.model.FamilyState;
import trainingpicker.model.TrainingMicrograph;
import trainingpicker.model.MicrographFamilyData;
import trainingpicker.model.MicrographFamilyState;
import trainingpicker.model.Particle;
import trainingpicker.model.TrainingParticle;
import trainingpicker.model.ParticlePicker;


public class UntiltedMicrographCanvas extends ImageCanvas implements
		MouseWheelListener {

	private ParticlePairPickerJFrame frame;
	private TiltedMicrographCanvas tiltedcanvas;
	private UntiltedParticle dragged;
	private ParticlePairPicker pppicker;
	private UntiltedMicrograph untiltedmicrograph;
	private ImageWindow iw;
	

	public UntiltedMicrographCanvas(ParticlePairPickerJFrame frame) {
		super(frame.getUntiltedMicrograph().getImage());
		this.untiltedmicrograph = frame.getUntiltedMicrograph();
		this.frame = frame;
		addMouseWheelListener(this);
		this.pppicker = frame.getParticlePairPicker();
		tiltedcanvas = new TiltedMicrographCanvas(frame);
		iw = new ImageWindow(imp, this);
		WindowUtils.centerScreen(0, iw);
	}
	
	public void updateMicrograph() {
		this.untiltedmicrograph = frame.getUntiltedMicrograph();
		imp = untiltedmicrograph.getImage();
		iw.setImage(imp);
		iw.setTitle(untiltedmicrograph.getName());
		tiltedcanvas.updateMicrograph();
		setImageUpdated();
		repaint();
	}
	
	
	public void mouseEntered(MouseEvent e) {
		if (frame.getTool() != Tool.PICKER) {
			super.mouseEntered(e);
			return;
		}
		setCursor(crosshairCursor);
	}
	
	public void mouseMoved(MouseEvent e) {
		if (frame.getTool() != Tool.PICKER) {
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

	public void mousePressed(MouseEvent e) {
		if (frame.getTool() != Tool.PICKER) {
			super.mousePressed(e);
			return;
		}
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (SwingUtilities.isRightMouseButton(e)) {
			setupScroll(x, y);
			tiltedcanvas.mousePressed(x, y);
		}
		UntiltedParticle p = untiltedmicrograph.getParticle(x, y, frame.getParticleSize());
		if (p != null) {
			if (SwingUtilities.isLeftMouseButton(e) && e.isControlDown()) {
				untiltedmicrograph.removeParticle(p);
				frame.updateMicrographsModel();
			} else if (SwingUtilities.isLeftMouseButton(e)) {
				dragged = p;
			}
		} else if (SwingUtilities.isLeftMouseButton(e)
				&& TrainingParticle.boxContainedOnImage(x, y, frame.getParticleSize(), imp)) {
			p = new UntiltedParticle(x, y, untiltedmicrograph);
			untiltedmicrograph.addParticle(p);
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
		dragged = null;
	}

	/**
	 * Updates particle position and repaints if onpick.
	 */
	@Override
	public void mouseDragged(MouseEvent e) {

		if (frame.getTool() != Tool.PICKER) {
			super.mouseDragged(e);
			return;
		}
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if (SwingUtilities.isRightMouseButton(e)) {
			scroll(e.getX(), e.getY());
			tiltedcanvas.mouseDragged(e.getX(), e.getY());
			return;
		}
		
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
		tiltedcanvas.mouseWheelMoved(x, y, rotation);

	}

	
	
	public TiltedMicrographCanvas getTiltedCanvas()
	{
		return tiltedcanvas;
	}
	
	public void paint(Graphics g) {
		super.paint(g);
		Graphics2D g2 = (Graphics2D) g;
		int x0 = (int) getSrcRect().getX();
		int y0 = (int) getSrcRect().getY();
		int index = 0;
		for (Particle p: untiltedmicrograph.getParticles())
		{
			drawShape(g2, p, x0, y0, index == untiltedmicrograph.getParticles().size() - 1);
			index ++;
		}
	}
	

	private void drawShape(Graphics2D g2, Particle p, int x0, int y0, boolean all) {
		int radius = frame.getParticleSize()/2;
		int x = (int) ((p.getX() - x0) * magnification);
		int y = (int) ((p.getY() - y0) * magnification);
		int distance = (int)(10 * magnification);

		if (frame.isShapeSelected(Shape.Rectangle) || all)
			g2.drawRect(x - radius, y - radius, radius * 2, radius * 2);
		if (frame.isShapeSelected(Shape.Circle) || all)
			g2.drawOval(x - radius, y - radius, radius * 2, radius * 2);
		if (frame.isShapeSelected(Shape.Center) || all) 
		{
			g2.drawLine(x, y - distance, x, y + distance);
			g2.drawLine(x + distance, y, x - distance, y);
		}
	}





}
