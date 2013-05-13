package xmipp.viewer.particlepicker.training;

import java.util.logging.Level;

import javax.swing.SwingUtilities;

import xmipp.utils.XmippDialog;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.SingleParticlePicker;
import xmipp.viewer.particlepicker.training.gui.SingleParticlePickerJFrame;
import xmipp.viewer.particlepicker.training.model.Mode;

class Main
{
	// 0 --> input metadata
	// 1 --> output dir
	// 2 --> mode

	// On Supervised
	// 3 --> number of threads for supervised mode
	// 4 --> fast mode for supervised mode
	// 5 --> incore for supervised mode

	// On Review
	// 3 -->external auto dir
	public static void main(String[] args)
	{

		final String[] myargs = args;

		SwingUtilities.invokeLater(new Runnable()
		{

			@Override
			public void run()
			{
				SingleParticlePickerJFrame tp = null;
				try
				{
					SingleParticlePicker ppicker = null;
					String selfile = myargs[0];
					String outputdir = myargs[1];
					//mode is the third argument 
					Mode mode = Mode.getMode(myargs[2]);
					
					if (mode == Mode.Manual || mode == Mode.ReadOnly)
						ppicker = new SingleParticlePicker(selfile, outputdir, mode);

					else if (mode == Mode.Supervised)
					{
						int index = 3;
						int threads = Integer.parseInt(myargs[index]);
						boolean fastmode = Boolean.parseBoolean(myargs[index + 1]);
						boolean incore = Boolean.parseBoolean(myargs[index + 2]);
						ppicker = new SingleParticlePicker(selfile, outputdir, threads, fastmode, incore);
					}

					else if (mode == Mode.Review)
					{
						String reviewfile = myargs[4];
						ppicker = new SingleParticlePicker(selfile, outputdir, reviewfile);
					}
			
					tp = new SingleParticlePickerJFrame(ppicker);
				}
				catch (Exception e)
				{
					System.out.println("Error catched on main");
					ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
					XmippDialog.showException(null, e);
					
					
				}
			}
		});

	}

}