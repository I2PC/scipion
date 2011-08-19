package model;

import java.util.logging.Level;
import javax.swing.SwingUtilities;

import gui.XmippParticlePickerJFrame;
import ij.ImagePlus;





class Main implements Runnable
{
	public static void main(String[] args)
	{
		SwingUtilities.invokeLater(new Main(args[0], args[1]));
	}
	private String outputdir;
	private String xmd;
	
	Main(String xmd, String outputdir)
	{
		this.xmd = xmd;
		this.outputdir = outputdir;
	}
	
	@Override
	public void run() {
		try {
			PPConfiguration.setMicrographsSelFile(xmd);
			PPConfiguration.setOutputDir(outputdir);
			XmippParticlePickerJFrame frame = new XmippParticlePickerJFrame();
		} catch (Exception e) {
			PPConfiguration.getLogger().log(Level.SEVERE, e.getMessage(), e);
		}
	}
}