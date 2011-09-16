package trainingpicker.model;

import ij.ImagePlus;
import ij.gui.Roi;
import ij.process.ImageProcessor;

import java.awt.Rectangle;
import java.util.List;

import javax.swing.ImageIcon;

import trainingpicker.gui.ParticleCanvas;

public class TrainingParticle extends Particle {
	
	protected Family family;
	protected ImagePlus img;
	protected ParticleCanvas canvas;
	protected double cost = 2;
	
	
	public TrainingParticle(int x, int y, Family family, TrainingMicrograph micrograph)
	{
		this(x, y, family, micrograph, 2);
	}
	
	public TrainingParticle(int x, int y, Family family, TrainingMicrograph micrograph, double cost)
	{
		super(x, y, micrograph);
		this.family = family;
		this.cost = cost;
	}
	
	public double getCost()
	{
		return cost;
	}


	public Family getFamily() {
		return family;
	}

	public void setFamily(Family family) {
		this.family = family;
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


	
	public ParticleCanvas getImageCanvas()
	{
		if(canvas == null)
		{
			canvas = new ParticleCanvas(getImage());
		}
		return canvas;
		
	}






	
}
