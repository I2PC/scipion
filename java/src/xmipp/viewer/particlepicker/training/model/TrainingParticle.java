package xmipp.viewer.particlepicker.training.model;

import ij.ImagePlus;
import ij.process.ImageProcessor;

import java.awt.Rectangle;

import javax.swing.ImageIcon;

import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.ImageGeneric;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.Family;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.PickerParticle;

public class TrainingParticle extends PickerParticle{
	
	protected Family family;
	protected ImagePlus img;
	protected double cost = 2;
	
	
	
	

	public TrainingParticle(int x, int y, Family family, Micrograph micrograph)
	{
		this(x, y, family, micrograph, 2);
	}
	
	public TrainingParticle(int x, int y, Family family, Micrograph micrograph, double cost)
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
		return super.contains(x2, y2, family.getSize());
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
			throw new IllegalArgumentException(XmippMessage.getOutOfBoundsMsg(String.format(" particle center: %s %s", x, y)));
		super.setPosition(x, y);
		
		img = null;
	}
	
	public ImageIcon getImageIcon()
	{
		
		ImageIcon icon = new ImageIcon(getImagePlus().getImage());
		return icon;
	}
	
	

	public ImageGeneric getImageGeneric() {
		try {
			return XmippImageConverter.convertToImageGeneric(getImagePlus());
		} catch (Exception e) {
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	
}
