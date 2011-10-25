package particlepicker;


import ij.ImagePlus;
import ij.gui.ImageCanvas;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

import javax.swing.SwingUtilities;
import particlepicker.training.model.TrainingParticle;

public abstract class ParticlePickerCanvas extends ImageCanvas implements MouseWheelListener
{


	public ParticlePickerCanvas(ImagePlus imp)
	{
		super(imp);
		// TODO Auto-generated constructor stub
	}

	public void moveTo(TrainingParticle p)
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
		if (!getSrcRect().contains(r))
		{
			setSourceRect(r);
			repaint();
		}
	}
	

	
	/**
	 * Adds particle or updates its position if onpick. If ondeletepick removes
	 * particle. Considers owner for selection to the first particle containing
	 * point. Sets dragged if onpick
	 */

	public void mousePressed(MouseEvent e)
	{
		if (getFrame().getTool() != Tool.PICKER)
		{
			super.mousePressed(e);
			return;
		}
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (SwingUtilities.isRightMouseButton(e))
		{
			setupScroll(x, y);
			return;
		}
		
	}
	
	public void mouseEntered(MouseEvent e)
	{
		if (getFrame().getTool() != Tool.PICKER)
		{
			super.mouseEntered(e);
			return;
		}
		setCursor(crosshairCursor);
	}

	public void mouseMoved(MouseEvent e)
	{
		if (getFrame().getTool() != Tool.PICKER)
		{
			super.mouseMoved(e);
			return;
		}
		setCursor(crosshairCursor);
	}
	
	
	public void mouseDragged(MouseEvent e)
	{

		if (getFrame().getTool() != Tool.PICKER)
		{
			super.mouseDragged(e);
			return;
		}
		if (SwingUtilities.isRightMouseButton(e))
		{
			scroll(e.getX(), e.getY());
			return;
		}
	}
	
	public void mouseReleased(MouseEvent e)
	{
		if (getFrame().getTool() != Tool.PICKER)
		{
			super.mouseReleased(e);
			return;
		}
		
	}
	
	public abstract void setActive(TrainingParticle p);
	
	public abstract ParticlePickerJFrame getFrame();
	
	
	protected void drawShape(Graphics2D g2, TrainingParticle p, boolean all)
	{
		int x0 = (int) getSrcRect().getX();
		int y0 = (int) getSrcRect().getY();
		int size = (int) (p.getFamily().getSize() * magnification);
		int radius = (int) (p.getFamily().getSize() / 2 * magnification);
		int x = (int) ((p.getX() - x0) * magnification);
		int y = (int) ((p.getY() - y0) * magnification);
		int distance = (int) (5 * magnification);

		if (getFrame().isShapeSelected(Shape.Rectangle) || all)
			g2.drawRect(x - radius, y - radius, size, size);
		if (getFrame().isShapeSelected(Shape.Circle) || all)
			g2.drawOval(x - radius, y - radius, size, size);
		if (getFrame().isShapeSelected(Shape.Center) || all)
		{
			g2.drawLine(x, y - distance, x, y + distance);
			g2.drawLine(x + distance, y, x - distance, y);
		}
	}
	
	@Override
	public void mouseWheelMoved(MouseWheelEvent e)
	{
		int x = e.getX();
		int y = e.getY();

		int rotation = e.getWheelRotation();
		if (rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);
		if (getMagnification() <= 1.0)
			imp.repaintWindow();

	}


}
