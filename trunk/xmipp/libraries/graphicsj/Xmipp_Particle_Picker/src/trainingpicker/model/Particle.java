package trainingpicker.model;

import ij.ImagePlus;
import ij.process.ImageProcessor;

import java.awt.Rectangle;
import java.util.List;

import javax.swing.ImageIcon;

import trainingpicker.gui.ParticleCanvas;

public class Particle {
	protected int x;
	protected int y;
	protected Micrograph micrograph;
	
	
	public Particle(int x, int y, Micrograph micrograph)
	{
		this.x = x;
		this.y = y;
		this.micrograph = micrograph;
	}
	
	
	public Micrograph getMicrograph() {
		return micrograph;
	}


	public void setMicrograph(TrainingMicrograph micrograph) {
		this.micrograph = micrograph;
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
	
	public boolean contains(int x2, int y2, int size )
	{
		int radius = size/2;
			if(x2 < x - radius || x2 > x + radius)
				return false;
			if(y2 < y - radius || y2 > y + radius)
				return false;
			return true;
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

	public void setPosition(int x, int y) {
		this.x = x;
		this.y = y;
		
	}
	
	public String toString()
	{
		return String.format("x = %s; y = %s", x, y); 
	}
	

	


}
