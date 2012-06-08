package xmipp.particlepicker.training.model;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;

import xmipp.particlepicker.Family;
import xmipp.utils.XmippMessage;
import xmipp.jni.Particle;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;

public class ReviewParticlePicker extends TrainingPicker
{

	private String reviewfile;
	private Family reviewfamily;

	public String getReviewFile()
	{
		return reviewfile;
	}

	public Family getReviewFamily()
	{
		return reviewfamily;
	}

	public ReviewParticlePicker(String selfile, String outputdir, String reviewfile)
	{
		super(selfile, outputdir, FamilyState.Review);
		if (!new File(reviewfile).exists())
			throw new IllegalArgumentException(XmippMessage.getNoSuchFieldValueMsg("review file", reviewfile));
		this.reviewfile = reviewfile;
		String[] parts = reviewfile.split(File.separator);
		String familyname = parts[parts.length - 1].split("_")[0];
		this.reviewfamily = getFamily(familyname);
		if (reviewfamily == null)
			throw new IllegalArgumentException(XmippMessage.getNoSuchFieldValueMsg("family", familyname));
		reviewfamily.setState(FamilyState.Review);
		loadMicrographs();
		families.clear();
		families.add(reviewfamily);
	}

	@Override
	public void persistMicrographs()
	{
		exportParticles(reviewfamily, reviewfile);
	}

	@Override
	public void loadMicrographs()
	{
		micrographs.clear();
		TrainingMicrograph micrograph;
		String ctf = null, filename;
		try
		{
			MetaData md = new MetaData(getMicrographsSelFile());
			boolean existsctf = md.containsLabel(MDLabel.MDL_PSD_ENHANCED);
			long[] ids = md.findObjects();
			for (long id : ids)
			{

				filename = md.getValueString(MDLabel.MDL_MICROGRAPH, id);
				if (existsctf)
					ctf = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				micrograph = new TrainingMicrograph(filename, ctf, families, getMode());
				micrographs.add(micrograph);
			}
			if (micrographs.size() == 0)
				throw new IllegalArgumentException(String.format("No micrographs specified on %s", getMicrographsSelFile()));
			importParticlesFromXmipp30Folder(reviewfamily, reviewfile);
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}
}
