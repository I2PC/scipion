package pairpicker.model;

import java.io.File;
import java.util.List;

import javax.swing.Icon;

import picker.model.Constants;
import picker.model.Family;
import picker.model.FamilyState;
import picker.model.MicrographFamilyData;
import picker.model.MicrographFamilyState;

public class MicrographPair {
	
	private String image;
	private String tiltedimage;
	private String name;
	private String tiltedname;
	
	public MicrographPair(String image, String tiltedimage) {
		
		this.image = image;
		this.tiltedimage = tiltedimage;
		if(!new File(image).exists())
			throw new IllegalArgumentException(Constants.getNoSuchFieldValueMsg("image", image));
		if(!new File(tiltedimage).exists())
			throw new IllegalArgumentException(Constants.getNoSuchFieldValueMsg("tilted image", tiltedimage));
		this.name = getName(image);
		this.tiltedname = getName(tiltedimage);

	}

	public void releaseImage() {
		// TODO Auto-generated method stub
		
	}
	
	public String getImage() {
		return image;
	}

	public String getTiltedimage() {
		return tiltedimage;
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

	public String getTiltedName() {
		// TODO Auto-generated method stub
		return tiltedname;
	}

}
