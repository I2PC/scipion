package gui;

import ij.ImagePlus;
import ij.gui.ImageCanvas;

import java.awt.BasicStroke;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Stroke;
import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.List;

import javax.swing.SwingUtilities;

import browser.windows.iPollImageWindow;

import model.Family;
import model.Micrograph;
import model.Particle;

public class PPCanvas extends ImageCanvas implements MouseWheelListener{

	private XmippParticlePickerJFrame frame;
	private Micrograph micrograph;
	private Particle dragged;
	final static BasicStroke dashedst = new BasicStroke(1.0f,
            BasicStroke.CAP_BUTT,
            BasicStroke.JOIN_MITER,
            10.0f, new float[]{10.0f}, 0.0f);
	final static BasicStroke continuousst = new BasicStroke();
	

	public PPCanvas(XmippParticlePickerJFrame frame, Micrograph micrograph) {
		super(micrograph.getImage());
		this.micrograph = micrograph;
		this.frame = frame;
		addMouseWheelListener(this);
		
		// TODO Auto-generated constructor stub
	}


	
	public Particle getParticle(int x, int y)
	{
		for(Particle p: micrograph.getParticles())
		{
			if (p.contains(x, y)) 
			return p;
		}
		return null;
	}
	
	/**
	 * Adds particle or updates its position if onpick.
	 * If ondeletepick removes particle. Considers owner for
	 * selection to the first particle containing point.
	 * Sets dragged if onpick
	 */
	public void mousePressed(MouseEvent e) {
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if(frame.getTool() != Tool.PICKER)
		{
			super.mousePressed(e);
			return;
		}
		if(SwingUtilities.isRightMouseButton(e))
		{
			setupScroll(x, y);
			return;
		}
		Particle p = getParticle(x, y);
		if (p != null)
		{
			if(SwingUtilities.isLeftMouseButton(e) && e.isControlDown())
			{
				micrograph.removeParticle(p);
				frame.updateMicrographsModel();
			}
			else if(SwingUtilities.isLeftMouseButton(e))
			{
				p.setPosition(x, y);
				dragged = p;
			}
		}
		else if (SwingUtilities.isLeftMouseButton(e)
				&& Particle.boxContainedOnImage(x, y, frame.getFamily().getSize(), imp)) {
			p = new Particle(x, y, frame.getFamily(), micrograph);
			micrograph.addParticle(p);
			dragged = p;
			frame.updateMicrographsModel();
		}
		frame.setChanged(true);
		repaint();
	}

	
	
	/**
	 * Updates particle position and repaints. 
	 * Sets dragged to null at the end
	 */
	public void mouseReleased(MouseEvent e) {
		if(frame.getTool() != Tool.PICKER)
		{
			super.mouseReleased(e);
			return;
		}
		
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if (dragged == null)//not onpick
			return;
		Particle p = null;
		p = dragged;
		dragged = null;
		if(!Particle.boxContainedOnImage(x, y, p.getFamily().getSize(), imp))
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
		if(frame.getTool() != Tool.PICKER)
		{
			super.mouseDragged(e);
			return;
		}
		if(SwingUtilities.isRightMouseButton(e))
		{
			scroll(e.getX(), e.getY());
			return;
		}
		if(dragged == null)
			return;
		
		if(!Particle.boxContainedOnImage(x, y, dragged.getFamily().getSize(), imp))
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
		
		for(Particle p: micrograph.getParticles())
		{
			g2.setColor(p.getFamily().getColor());
			if(p.isAuto())
				g2.setStroke(dashedst);
			else
				g2.setStroke(continuousst);
			radius = (int)(p.getFamily().getSize()/ 2 * magnification);
			count++;
			x = (int) ((p.getX() - x0) * magnification);
			y = (int) ((p.getY() - y0) * magnification);
			g2.drawString(Integer.toString(count), x, y - radius-2);
			if(frame.getShape() == Shape.RECTANGLE || frame.getShape() == Shape.BOTH)
				g2.drawRect(x - radius , y - radius, radius * 2, radius * 2);
			if (frame.getShape() == Shape.CIRCLE || frame.getShape() == Shape.BOTH)
				g2.drawOval(x - radius , y - radius, radius * 2, radius * 2);
			
		}
	}
	
	public void setMicrograph(Micrograph micrograph)
	{
		this.micrograph = micrograph;
		imp = micrograph.getImage();
		setImageUpdated();
		repaint();
	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent e) {
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		
		int rotation = e.getWheelRotation();
		if(rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);
		
	}
	
	
	

}
