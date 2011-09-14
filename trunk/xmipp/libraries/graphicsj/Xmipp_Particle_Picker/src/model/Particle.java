package model;

import gui.ParticleCanvas;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.process.ImageProcessor;

import java.awt.Rectangle;
import java.util.List;

import javax.swing.ImageIcon;

public class Particle implements Comparable<Particle>{
	
	protected int x;
	protected int y;
	protected Family family;
	protected Micrograph micrograph;
	protected ImagePlus img;
	protected ParticleCanvas canvas;
	protected double cost = 2;
	
	
	public Particle(int x, int y, Family family, Micrograph micrograph)
	{
		this(x, y, family, micrograph, 2);
	}
	
	public Particle(int x, int y, Family family, Micrograph micrograph, double cost)
	{
		this.x = x;
		this.y = y;
		this.family = family;
		this.micrograph = micrograph;
		this.cost = cost;
	}
	
	public double getCost()
	{
		return cost;
	}
	
	
	public Micrograph getMicrograph() {
		return micrograph;
	}


	public void setMicrograph(Micrograph micrograph) {
		this.micrograph = micrograph;
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
	
	public boolean contains(int x2, int y2 )
	{
		int radius = family.getSize()/2;
			if(x2 < x - radius || x2 > x + radius)
				return false;
			if(y2 < y - radius || y2 > y + radius)
				return false;
			return true;
	}



	
	public ImagePlus getImage()
	{
		if(img == null)
		{
			int size = family.getSize();
			ImagePlus mimage = micrograph.getImage();
			int radius = size/2;
			Rectangle r = new Rectangle(x - radius , y - radius, radius * 2, radius * 2);
			mimage.setRoi(r);
			ImageProcessor processor = mimage.getProcessor().crop();
			img = new ImagePlus("", processor);
		}
		return img;
	}
	
	public ImageIcon getImageIcon()
	{
		
		ImageIcon icon = new ImageIcon(img.getImage());
		return icon;
	}

	public static boolean boxContainedOnImage(int x, int y, int size, ImagePlus img) {
		if(img == null)
			return false;
		int width = img.getWidth();
		int height = img.getHeight();
		int radius = size/2;
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

	
	public ParticleCanvas getImageCanvas()
	{
		if(canvas == null)
		{
			canvas = new ParticleCanvas(getImage());
		}
		return canvas;
		
	}


	public void setPosition(int x, int y) {
		this.x = x;
		this.y = y;
		
	}
	
	public String toString()
	{
		return String.format("x = %s; y = %s", x, y); 
	}




	@Override
	public int compareTo(Particle p) {
		if(p.x > x)
			return 1;
		if(p.x == x)
		{
			if(p.y > y)
				return 1;
			if(p.y == y)
				return 0;
			return -1;
		}
		return -1;
	}


	
}
