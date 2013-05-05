package xmipp.viewer.particlepicker.training.model;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;

import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.Family;
import xmipp.jni.Particle;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;

public class ReviewParticlePicker extends TrainingPicker
{

	private String reviewfile;

	public String getReviewFile()
	{
		return reviewfile;
	}


	
	public ReviewParticlePicker(String selfile, String outputdir, String fname, String reviewfile)
	{
		super(selfile, outputdir, fname, Mode.Review);
		if (!new File(reviewfile).exists())
			throw new IllegalArgumentException(XmippMessage.getNoSuchFieldValueMsg("review file", reviewfile));
		this.reviewfile = reviewfile;
		importAllParticles(reviewfile);
		
	}

	public ReviewParticlePicker(String selfile, String outputdir, String reviewfile)
	{
		this(selfile, outputdir, getFamilyName(reviewfile), reviewfile);
	
	}
	
	public static String getFamilyName(String reviewfile)
	{
		String[] parts = reviewfile.split(File.separator);
		String familyname = parts[parts.length - 1].split("_")[0];
		return familyname;
	}

	@Override
	public void saveMicrographs()
	{
		exportParticles(reviewfile);
	}
	
	

	@Override
	public void loadMicrographs()
	{
		
		try
		{
			loadEmptyMicrographs();
			importAllParticles(reviewfile);
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}



	
}
