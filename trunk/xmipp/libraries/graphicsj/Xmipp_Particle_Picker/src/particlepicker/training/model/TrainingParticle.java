package particlepicker.training.model;

import ij.ImagePlus;
import ij.process.ImageProcessor;

import java.awt.Rectangle;

import javax.swing.ImageIcon;

import particlepicker.Constants;
import particlepicker.Family;
import particlepicker.Micrograph;
import particlepicker.ParticleCanvas;
import particlepicker.ParticlePickerJFrame;
import particlepicker.training.gui.TrainingPickerJFrame;
import xmipp.Particle;

public class TrainingParticle extends Particle{
	
	protected Family family;
	protected ImagePlus img;
	protected double cost = 2;
	protected Micrograph micrograph;
	private ParticleCanvas canvas;
	

	public TrainingParticle(int x, int y, Family family, Micrograph micrograph)
	{
		this(x, y, family, micrograph, 2);
	}
	
	public TrainingParticle(int x, int y, Family family, Micrograph micrograph, double cost)
	{
		super(x, y);
		this.micrograph = micrograph;
		this.family = family;
		this.cost = cost;
	}
	
	public Micrograph getMicrograph() {
		return micrograph;
	}


	public void setMicrograph(Micrograph micrograph) {
		this.micrograph = micrograph;
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
	
	public int getX0()
	{
		return getX() - family.getRadius();
	}
	
	public int getY0()
	{
		return getY() - family.getRadius();
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

	
	public ImagePlus getImagePlus()
	{
		if(img == null)
		{
			int size = family.getSize();
			ImagePlus mimage = micrograph.getImagePlus();
			int radius = size/2;
			Rectangle r = new Rectangle(x - radius , y - radius, radius * 2, radius * 2);
			mimage.setRoi(r);
			ImageProcessor processor = mimage.getProcessor().crop();
			img = new ImagePlus("", processor);
		}
		return img;
	}
	
	@Override
	public void setPosition(int x, int y)
	{
		int radius = family.getSize()/2;
		if(x - radius < 0 || y - radius < 0 || x + radius > micrograph.getImagePlus().getWidth() || y + radius > micrograph.getImagePlus().getHeight())
			throw new IllegalArgumentException(Constants.getOutOfBoundsMsg(String.format(" particle center: %s %s", x, y)));
		super.setPosition(x, y);
		
		img = null;
	}
	
	public ImageIcon getImageIcon()
	{
		
		ImageIcon icon = new ImageIcon(getImagePlus().getImage());
		return icon;
	}
	
	public ParticleCanvas getParticleCanvas(ParticlePickerJFrame frame)
	{
		if(canvas == null)
			canvas = new ParticleCanvas(this, frame);
		return canvas; 
	}
	
	public void resetParticleCanvas()
	{
		canvas = null;
	}

	
}
