package model;

import java.io.File;
import java.io.IOException;
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

public class PPConfiguration {
	
	private static Logger logger;
	private static String outputdir = ".";
	private static String mgxmd = "micrographs.sel";
	private static String rundir;
	
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
	
	
	public static String getOutputDir()
	{
		return outputdir;
	}
	
	public static String getFamiliesXMD()
	{
		return "families.xmd";
	}
	
	public static String getMicrographsXMD()
	{
		return getMicrographPath(mgxmd);
	}
	
	public static String getMicrographPath(String rpath)
	{
		return rundir + File.separator + rpath;
	}
	
	public static void setMicrographsXMD(String xmd)
	{
		mgxmd = xmd;
	}
	
	public static void setRunDir(String dir)
	{
		rundir = dir;
	}
	
	public static String getRunDir()
	{
		return rundir;
	}
	

}
