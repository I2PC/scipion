package tiltpairpicker.model;

import ij.ImagePlus;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import trainingpicker.model.Constants;
import trainingpicker.model.Particle;

public class UntiltedMicrograph {
	
	private String imagefile;
	private TiltedMicrograph tiltedmicrograph;
	private String name;
	private ImagePlus image;
	private List<Particle> particles;
	
	public UntiltedMicrograph(String image, TiltedMicrograph tiltedmicrograph) {
		
		this.imagefile = image;
		this.tiltedmicrograph = tiltedmicrograph;
		if(!new File(image).exists())
			throw new IllegalArgumentException(Constants.getNoSuchFieldValueMsg("image", image));
		this.name = getName(image);
		particles = new ArrayList<Particle>();
	}

	
	public String getImageFile() {
		return imagefile;
	}

	public String getName() {
		return name;
	}
	
	public static String getName(String file)
	{
		String[] tokens = file.split(File.separator);
		return  tokens[tokens.length - 1];
	}
	
	public void setName(String name) {
		this.name = name;
	}

	public boolean isEmpty() {
		// TODO Auto-generated method stub
		return false;
	}

	public ImagePlus getImage()
	{
		if(image == null)
			image = new ImagePlus(imagefile);
		return image;
	}
	
	public TiltedMicrograph getTiltedMicrograph()
	{
		return tiltedmicrograph;
	}
	
	public void releaseImage()
	{
		image = null;
	}

}
