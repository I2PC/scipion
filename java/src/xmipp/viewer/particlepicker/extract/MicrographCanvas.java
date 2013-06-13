package xmipp.viewer.particlepicker.extract;

import ij.gui.ImageWindow;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.util.List;
import xmipp.ij.commons.XmippImageCanvas;
import xmipp.jni.Particle;
import xmipp.viewer.particlepicker.Shape;

public class MicrographCanvas extends XmippImageCanvas
{

	
	private List<Particle> particles;
	private int size;
	private MicrographData micdata;

	public MicrographCanvas(MicrographData loader, int size)
	{
		super(loader.getMicrograph().getImagePlus());
		this.micdata = loader;
		this.particles = loader.getParticles();
		this.size = size;
		// TODO Auto-generated constructor stub
	}
	
	protected ImageWindow iw;
	
	public void display()
	{
		if (iw != null && iw.isVisible())
		{
			iw.setImage(getImage());
			iw.updateImage(getImage());
		}
		else
			this.iw = new ImageWindow(getImage(), this);// if you dont provide
														// iw, I init mine
		// iw.maximize();
		iw.pack();
	}

	public void paint(Graphics g)
	{
		Graphics offgc;//Off screen graphic
		Image offscreen = null;//Off screen image
		Dimension d = getSize();
		// create the offscreen buffer and associated Graphics
		offscreen = createImage(d.width, d.height);
		offgc = offscreen.getGraphics();
		super.paint(offgc);//super paint in off screen
		//my paint in offscreen
		Graphics2D g2 = (Graphics2D) offgc;
		doCustomPaint(g2);
		//drawing offscreen image
		g.drawImage(offscreen, 0, 0, this);
	}

	protected void drawShape(Graphics2D g2, Shape shape, int x, int y, int radius)
	{
		int size = 2 * radius;
		if (shape == Shape.Rectangle)
			g2.drawRect(x - radius, y - radius, size, size);
		else if (shape == Shape.Circle)
			g2.drawOval(x - radius, y - radius, size, size);
		else if (shape == Shape.Center)
		{
			int distance = (int) (radius / 5. * magnification);
			g2.drawLine(x, y - distance, x, y + distance);
			g2.drawLine(x + distance, y, x - distance, y);
		}

	}

	protected void doCustomPaint(Graphics2D g2)
	{
		g2.setColor(Color.GREEN);
		for (Particle p : particles)
		{
			drawShape(g2, Shape.Rectangle, getXOnImage(p.getX()), getYOnImage(p.getY()), size);
			drawShape(g2, Shape.Center, getXOnImage(p.getX()), getYOnImage(p.getY()), size);
		}
	}

}
