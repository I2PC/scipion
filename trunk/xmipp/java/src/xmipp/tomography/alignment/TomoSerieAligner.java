package xmipp.tomography.alignment;

import ij.ImageStack;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;

import xmipp.utils.XmippMessage;

public class TomoSerieAligner {
	
	public static void main(String[] args)
	{
		TomoSerieAligner aligner = new TomoSerieAligner(args[0], args[1], 50, 45);
		aligner.alignTomographies();
	}
	
	public static Logger getLogger()
	{
		try
		{
			if (logger == null)
			{
				FileHandler fh = new FileHandler("TomoSerieAligner.log", true);
				fh.setFormatter(new SimpleFormatter());
				logger = Logger.getLogger("TomoSerieAlignerLogger");
				logger.addHandler(fh);
			}
			return logger;
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	
	private static Logger logger;
	private String tomoseriefile;
	private String outputdir;
	private List<Tomography> tomographies;
	private int dmin;
	private int dmax;
	
	public TomoSerieAligner(String tomoseriefile, String outputdir, int dmin, int dmax)
	{
		this.tomoseriefile = tomoseriefile;
		this.outputdir = outputdir;
		this.dmin = dmin;
		this.dmax = dmax;
		loadTomographies();
	}
	
	public int getDmin() {
		return dmin;
	}



	public int getDmax() {
		return dmax;
	}



	private void loadTomographies()
	{
		tomographies = new ArrayList<Tomography>();
		String file = tomoseriefile;
		if (!new File(file).exists())
			throw new IllegalArgumentException(XmippMessage.getNoSuchFieldValueMsg("tomoseriefile", tomoseriefile));
		Tomography tomography = null, previous;
		String tomofile;
		int tiltangle;
		try
		{
			MetaData md = new MetaData(file);
			long[] ids = md.findObjects();
			for (long id : ids)
			{
				tomofile = md.getValueString(MDLabel.MDL_IMAGE, id);
//				System.out.println(tomofile);
				tiltangle = md.getValueInt(MDLabel.MDL_ANGLE_TILT, id);
				previous = tomography;
				tomography = new Tomography(tomofile, tiltangle, previous);
				if(previous != null)
					previous.setNext(tomography);
				tomographies.add(tomography);
			}
			if (tomographies.isEmpty())
				throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("tomographies"));
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}
	
	public List<Tomography> getTomographies()
	{
		return tomographies;
	}
	
	public ImageStack alignTomographies()
	{
		System.out.println(new Date());
		for(Tomography t: tomographies)
		{
			if(t.getPrevious() != null)
			{
				t.computeAffineTransform();
			}
		}
		System.out.println(new Date());
		return null;
	}
	

}
