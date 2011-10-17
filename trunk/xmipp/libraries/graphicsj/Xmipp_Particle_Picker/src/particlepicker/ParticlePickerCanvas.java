package particlepicker;

import ij.ImagePlus;
import ij.gui.ImageCanvas;

import java.awt.Graphics2D;
import java.awt.Rectangle;

import particlepicker.training.model.TrainingParticle;

public abstract class ParticlePickerCanvas extends ImageCanvas
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


}
