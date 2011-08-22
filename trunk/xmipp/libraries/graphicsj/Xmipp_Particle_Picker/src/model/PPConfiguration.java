package model;

import java.io.File;
import java.io.IOException;
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

public class PPConfiguration {
	
	private static Logger logger;
	private static String outputdir = ".";
	private static String mgselfile = "micrographs.sel";
	private static String rundir = ".";
	private static int threads;
	private static boolean fastMode;
	private static boolean incore;
	private static boolean auto;
	private static int minparticles = 100;
	
	public static Logger getLogger()
	{
		try {
			if(logger == null)
			{
				FileHandler fh = new FileHandler("PPicker.log", true);
				fh.setFormatter(new SimpleFormatter());
				logger = Logger.getLogger("PPickerLogger");
				logger.addHandler(fh);
			}
			return logger;
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	public static boolean getIsAuto()
	{
		return auto;
	}
	
	public static int getMinParticles()
	{
		return minparticles;
	}
	
	public static String getOutputDir()
	{
		return outputdir;
	}
	
	
	public static String getMicrographPath(String rpath)
	{
		return rundir + File.separator + rpath;
	}
	
	public static void setMicrographsSelFile(String file)
	{
		mgselfile = file;
	}
	
	public static String getMicrographsSelFile()
	{
		return mgselfile;
	}
	
	public static void setOutputDir(String dir)
	{
		if(!new File(dir).exists())
			throw new IllegalArgumentException("Output dir " + dir + " does not exist");
		outputdir = dir;
	}
	
	public static String getOutputPath(String file)
	{
		return outputdir + File.separator + file;
	}


	public static void setThreads(int threads) {
		PPConfiguration.threads = threads;
		
	}


	public static void setFastMode(boolean fastMode) {
		PPConfiguration.fastMode = fastMode;
	}
	
	public static void setIncore(boolean incore) {
		PPConfiguration.incore = incore;
	}

	public static boolean isFastMode() {
		return fastMode;
	}

	public static boolean isIncore() {
		return incore;
	}
	
	public static boolean getFastMode()
	{
		return fastMode;
	}

	public static void setIsAuto(boolean automatic) {
		PPConfiguration.auto = automatic;
		
	}
	
	public static int getThreads()
	{
		return threads;
	}

}
