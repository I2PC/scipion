package model;

import java.util.logging.Level;
import javax.swing.SwingUtilities;

import gui.XmippParticlePickerJFrame;
import ij.ImagePlus;





class MyRunnable implements Runnable
{
	public static void main(String[] args)
	{
		SwingUtilities.invokeLater(new MyRunnable(args[0], args[1]));
	}
	private String rundir;
	private String xmd;
	
	MyRunnable(String rundir, String xmd)
	{
		this.rundir = rundir;
		this.xmd = xmd;
	}
	
	@Override
	public void run() {
		try {
			PPConfiguration.setRunDir(rundir);
			PPConfiguration.setMicrographsXMD(xmd);
			XmippParticlePickerJFrame frame = new XmippParticlePickerJFrame();
		} catch (Exception e) {
			PPConfiguration.getLogger().log(Level.SEVERE, e.getMessage(), e);
		}
	}
}