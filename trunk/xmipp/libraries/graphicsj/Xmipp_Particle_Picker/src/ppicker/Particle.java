package ppicker;

import ij.ImagePlus;
import ij.gui.Roi;
import ij.process.ImageProcessor;

import java.awt.Rectangle;
import java.util.List;

import javax.swing.ImageIcon;

public class Particle {
	
	private int x;
	private int y;
	private boolean dragged;
	private Family family;

	
	
	public Particle(int x, int y, Family family)
	{
		this.x = x;
		this.y = y;
		this.family = family;
	}


	public Family getFamily() {
		return family;
	}

	public void setFamily(Family family) {
		this.family = family;
	}

	public int getX() {
		return x;
	}

	public void setX(int x) {
		this.x = x;
	}

	public int getY() {
		return y;
	}

	public void setY(int y) {
		this.y = y;
	}
	
	public boolean contains(int radius, int x2, int y2 )
	{
			if(x2 < x - radius || x2 > x + radius)
				return false;
			if(y2 < y - radius || y2 > y + radius)
				return false;
			return true;
	}

	public boolean isDragged() {
		return dragged;
	}

	public void setDragged(boolean dragged) {
		this.dragged = dragged;
	}
	
	public ImagePlus getImage(ImagePlus container, int radius)
	{
		Rectangle r = new Rectangle(x - radius , y - radius, radius * 2, radius * 2);
		container.setRoi(r);
		ImageProcessor processor = container.getProcessor().crop();
		return new ImagePlus("", processor);
	}
	
	public ImageIcon getImageIcon(ImagePlus container, int radius)
	{
		ImagePlus img = getImage(container, radius);
		ImageIcon icon = new ImageIcon(img.getImage());
		return icon;
	}

	public static boolean contained(int x, int y, int radius, ImagePlus img) {
		if(img == null)
			return false;
		int width = img.getWidth();
		int height = img.getHeight();
		if(x - radius < 0)
			return false;
		if(x + radius > width)
			return false;
		if(y - radius < 0)
			return false;
		if(y + radius > height)
			return false;
		return true;
			
	}


	public static int getDefaultRadius() {
		// TODO Auto-generated method stub
		return 50;
	}

}
