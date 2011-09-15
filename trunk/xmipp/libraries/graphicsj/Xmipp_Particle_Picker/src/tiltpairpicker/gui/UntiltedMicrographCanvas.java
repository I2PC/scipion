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

import picker.gui.WindowUtils;
import picker.model.AutomaticParticle;
import picker.model.FamilyState;
import picker.model.Micrograph;
import picker.model.MicrographFamilyData;
import picker.model.MicrographFamilyState;
import picker.model.Particle;
import picker.model.ParticlePicker;
import tiltpairpicker.model.ParticlePairPicker;
import tiltpairpicker.model.UntiltedMicrograph;


public class UntiltedMicrographCanvas extends ImageCanvas implements
		MouseWheelListener {

	private ParticlePairPickerJFrame frame;
	private TiltedMicrographCanvas tiltedcanvas;
	private Particle dragged;
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

	public void paint(Graphics g) {
		super.paint(g);
		Graphics2D g2 = (Graphics2D) g;
		
	}
	
	

	private void drawShape(Graphics2D g2, Particle p, int x0, int y0, int radius, boolean all) {
		
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
