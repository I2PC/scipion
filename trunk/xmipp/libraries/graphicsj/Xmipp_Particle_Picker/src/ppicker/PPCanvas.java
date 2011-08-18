package ppicker;

import ij.ImagePlus;
import ij.gui.ImageCanvas;

import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.List;

public class PPCanvas extends ImageCanvas implements MouseWheelListener{

	private XmippParticlePickerJFrame frame;
	private Particle dragged;

	public PPCanvas(ImagePlus imp, XmippParticlePickerJFrame frame) {
		super(imp);
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
		List<Family> families = frame.getPPData().getFamilies();
		Family family;
		for (int i = 0; i < families.size(); i++) 
		{
			family = families.get(i);
			for(Particle p2: family.getParticles())
			if (p2.contains(family.getRadius(), x, y)) {
				p = p2;
				if((e.getModifiers() & InputEvent.BUTTON1_MASK) != 0)
				{
					p.setX(x);
					p.setY(y);
					dragged = p;
				}
				else if((e.getModifiers() & InputEvent.BUTTON3_MASK) != 0)
				{
					p.getFamily().removeParticle(p);
					break;
				}
			}
		}
		if (((e.getModifiers() & InputEvent.BUTTON1_MASK) != 0) 
				&& p == null && Particle.contained(x, y, frame.getFamily().getRadius(), imp)) {
			p = new Particle(x, y, frame.getFamily());
			frame.getFamily().addParticle(p);
			dragged = p;
		}
		
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
		if(!Particle.contained(x, y, p.getFamily().getRadius(), imp))
			return;
		p.setX(x);
		p.setY(y);
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
		if(!Particle.contained(x, y, dragged.getFamily().getRadius(), imp))
			return;
		dragged.setX(x);
		dragged.setY(y);
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
		List<Family> families = frame.getPPData().getFamilies();
		for (Family family : families)
		{
			g2.setColor(family.getColor());
			radius = (int)(family.getRadius() * magnification);
			for(Particle p: family.getParticles())
			{
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
