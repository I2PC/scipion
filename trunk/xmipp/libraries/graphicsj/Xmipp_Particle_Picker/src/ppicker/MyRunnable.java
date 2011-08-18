package ppicker;

import java.util.logging.Level;
import javax.swing.SwingUtilities;
import ij.ImagePlus;





class MyRunnable implements Runnable
{
	public static void main(String[] args)
	{
		SwingUtilities.invokeLater(new MyRunnable(args[0]));
	}
	private String arg;
	MyRunnable(String arg)
	{
		this.arg = arg;
	}
	
	@Override
	public void run() {
		try {
			ImagePlus img  = new ImagePlus(arg);
			XmippParticlePickerJFrame frame = new XmippParticlePickerJFrame(img);
		} catch (Exception e) {
			PPConfiguration.getLogger().log(Level.SEVERE, e.getMessage(), e);
		}
	}
}