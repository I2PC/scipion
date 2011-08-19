package Xmipp;

import java.util.logging.Level;
import javax.swing.SwingUtilities;

import model.PPConfiguration;

import gui.XmippParticlePickerJFrame;
import ij.ImagePlus;





class Main implements Runnable
{
	//0 --> input metadata
	//1 --> output dir
	//2 --> number of threads for automatic picking
	//3 --> fast mode for automatic picking
	//4 --> incore for automatic picking
	public static void main(String[] args)
	{
		String mgselfile = args[0];
		String outputdir = args[1];
		
		if(args.length == 2)
			SwingUtilities.invokeLater(new Main(mgselfile, outputdir));
		else
		{
			int threads = Integer.parseInt(args[2]);
			boolean fastMode = Boolean.parseBoolean(args[3]);
			boolean incore = Boolean.parseBoolean(args[4]);
			SwingUtilities.invokeLater(new Main(mgselfile, outputdir, threads, fastMode, incore));
		}
			
	}
	
	
	Main(String mgselfile, String outputdir)
	{
		PPConfiguration.setMicrographsSelFile(mgselfile);
		PPConfiguration.setOutputDir(outputdir);
		PPConfiguration.setIsAuto(false);
	}
	
	Main(String mgselfile, String outputdir, int threads, boolean fastMode, boolean incore)
	{
		
		PPConfiguration.setMicrographsSelFile(mgselfile);
		PPConfiguration.setOutputDir(outputdir);
		PPConfiguration.setThreads(threads);
		PPConfiguration.setFastMode(fastMode);
		PPConfiguration.setIncore(incore);
		PPConfiguration.setIsAuto(true);
	}
	
	@Override
	public void run() {
		try {
			
			XmippParticlePickerJFrame frame = new XmippParticlePickerJFrame();
			
		} catch (Exception e) {
			PPConfiguration.getLogger().log(Level.SEVERE, e.getMessage(), e);
		}
	}
}