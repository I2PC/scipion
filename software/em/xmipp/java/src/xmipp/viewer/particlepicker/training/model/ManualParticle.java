package xmipp.viewer.particlepicker.training.model;

import ij.ImagePlus;
import ij.gui.Roi;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.ImageIcon;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.ImageGeneric;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.PickerParticle;

public class ManualParticle extends PickerParticle{
	
	protected ParticlePicker picker;
    protected ImageGeneric ig;

	protected double[] lastalign;

	public double[] getLastalign()
	{
		return lastalign;
	}

	public void setLastalign(double[] lastalign)
	{
		this.lastalign = lastalign;
	}

	public ManualParticle(int x, int y, ParticlePicker picker, Micrograph micrograph)
	{
		this(x, y, picker, micrograph, 2);
	}
	
	public ManualParticle(int x, int y, ParticlePicker picker, Micrograph micrograph, double cost)
	{
		super(x, y, micrograph, cost);
		this.picker = picker;
	}
	

	public ParticlePicker getParticlePicker()
	{
		return picker;
	}

	
	public int getX0()
	{
		return getX() - picker.getRadius();
	}
	
	public int getY0()
	{
		return getY() - picker.getRadius();
	}
	

	
	public boolean contains(int x2, int y2 )
	{
		return super.contains(x2, y2, picker.getSize());
	}

	
	protected void loadImagePlus()
	{
		if(ig == null)
		{
                    try {
                        int size = picker.getSize();
                        ImagePlus mimage = micrograph.getImagePlus();
                        int radius = size/2;
                        Rectangle r = new Rectangle(x - radius , y - radius, size, size);
                        Roi roi = mimage.getRoi();
                        mimage.setRoi(r);
                        ImageProcessor processor = mimage.getProcessor().crop();
                        ImagePlus img = new ImagePlus("", processor);
                        mimage.setRoi(roi);
                        ig = XmippImageConverter.convertToImageGeneric(img);
                    } catch (Exception ex) {
                        Logger.getLogger(ManualParticle.class.getName()).log(Level.SEVERE, null, ex);
                        throw new IllegalArgumentException(ex.getMessage());
                    }
		}
	}
	
	@Override
	public void setPosition(int x, int y)
	{
		int radius = picker.getSize()/2;
		if(!getMicrograph().fits(x, y, picker.getSize()))
			throw new IllegalArgumentException(XmippMessage.getOutOfBoundsMsg(String.format("Particle centered at %s, %s with size %s", x, y, picker.getSize())));
		super.setPosition(x, y);
		ig = null;
	}
	
	
	
	public void resetImage()
	{
//                if(ig != null)
//                    ig.destroy();
                lastalign = null;
		ig = null;
	}
	

	public ImageGeneric getImageGeneric() {
                loadImagePlus();//to ensure there is an image generic setted
                return ig;
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

	public int getSize()
	{
		return picker.getSize();
	}


	
}
