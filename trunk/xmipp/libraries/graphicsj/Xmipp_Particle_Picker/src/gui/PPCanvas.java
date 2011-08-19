package gui;

import ij.ImagePlus;
import ij.gui.ImageCanvas;

import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.List;

import model.Family;
import model.Micrograph;
import model.Particle;

public class PPCanvas extends ImageCanvas implements MouseWheelListener{

	private XmippParticlePickerJFrame frame;
	private Micrograph micrograph;
	private Particle dragged;

	public PPCanvas(XmippParticlePickerJFrame frame, Micrograph micrograph) {
		super(micrograph.getImage());
		this.micrograph = micrograph;
		this.frame = frame;
		addMouseWheelListener(this);
		
		// TODO Auto-generated constructor stub
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
		
		Particle p = null;
		Family family;
		for(Particle p2: micrograph.getParticles())
		{
			family = p2.getFamily();
			if (p2.contains(family.getSize(), x, y)) {
				p = p2;
				if((e.getModifiers() & InputEvent.BUTTON1_MASK) != 0)
				{
					p.setX(x);
					p.setY(y);
					dragged = p;
				}
				else if((e.getModifiers() & InputEvent.BUTTON3_MASK) != 0)
				{
					micrograph.removeParticle(p);
					frame.updateMicrographsModel();
					break;
				}
			}
		}
		if (((e.getModifiers() & InputEvent.BUTTON1_MASK) != 0) 
				&& p == null && Particle.boxContainedOnImage(x, y, frame.getFamily().getSize(), imp)) {
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
		p.setX(x);
		p.setY(y);
		frame.setChanged(true);
		repaint();
	}
	
		
	/**
	 * Updates particle position and repaints if onpick.
	 */
	@Override
	public void mouseDragged(MouseEvent e) {
		if(frame.getTool() != Tool.PICKER)
		{
			super.mouseDragged(e);
			return;
		}
		if(dragged == null)
			return;
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if(!Particle.boxContainedOnImage(x, y, dragged.getFamily().getSize(), imp))
			return;
		dragged.setX(x);
		dragged.setY(y);
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
