package xmipp.viewer.particlepicker.training.model;

import ij.ImagePlus;
import ij.gui.Roi;
import ij.process.ImageProcessor;

import java.awt.Rectangle;

import javax.swing.ImageIcon;

import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.Family;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticleCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.training.gui.TrainingPickerJFrame;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.ImageGeneric;
import xmipp.jni.Particle;

public class TrainingParticle extends Particle{
	
	protected Family family;
	protected ImagePlus img;
	protected double cost = 2;
	protected Micrograph micrograph;
	private ParticleCanvas canvas;
	protected double[] lastalign;
	
	
	

	public double[] getLastalign()
	{
		return lastalign;
	}

	public void setLastalign(double[] lastalign)
	{
		this.lastalign = lastalign;
	}
	
	

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
			Rectangle r = new Rectangle(x - radius , y - radius, size, size);
			Roi roi = mimage.getRoi();
			mimage.setRoi(r);
			ImageProcessor processor = mimage.getProcessor().crop();
			img = new ImagePlus("", processor);
			mimage.setRoi(roi);
		}
		return img;
	}
	
	@Override
	public void setPosition(int x, int y)
	{
		int radius = family.getSize()/2;
		if(!getMicrograph().fits(x, y, getFamily().getSize()))
			throw new IllegalArgumentException(XmippMessage.getOutOfBoundsMsg(String.format("Particle centered at %s, %s with size %s", x, y, getFamily().getSize())));
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

	public ImageGeneric getImageGeneric() {
		try {
			return XmippImageConverter.convertToImageGeneric(getImagePlus());
		} catch (Exception e) {
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void resetImagePlus()
	{
		img = null;
		
	}
	public double getTemplateRotation()
	{
		if(lastalign == null)
			return -1;
		return lastalign[1];
	}
	
	public double getTemplateTilt()
	{
		if(lastalign == null)
			return -1;
		return lastalign[2];
	}
	
	public double getTemplatePsi()
	{
		if(lastalign == null)
			return -1;
		return lastalign[3];
	}
	
	public int getTemplateIndex()
	{
		if(lastalign == null)
			return -1;
		return (int)lastalign[0];
	}

	
}
