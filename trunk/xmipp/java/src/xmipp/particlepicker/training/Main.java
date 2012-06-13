package xmipp.particlepicker.training;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import xmipp.particlepicker.training.gui.TrainingPickerJFrame;
import xmipp.particlepicker.training.model.FamilyState;
import xmipp.particlepicker.training.model.ManualParticlePicker;
import xmipp.particlepicker.training.model.ReviewParticlePicker;
import xmipp.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.particlepicker.training.model.TrainingPicker;

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
				TrainingPickerJFrame tp = null;
				try
				{
					TrainingPicker ppicker = null;
					String selfile = myargs[0];
					String outputdir = myargs[1];
					FamilyState mode;
					String fname = null;
					if(myargs.length == 3)
						mode = FamilyState.getFamilyState(myargs[2]);
					else
					{
						fname = myargs[2];
						mode = FamilyState.getFamilyState(myargs[3]);
					}

					if (mode == FamilyState.Manual)
						ppicker = (fname == null)? new ManualParticlePicker(selfile, outputdir, mode): new ManualParticlePicker(selfile, outputdir, fname, mode);

					else if (mode == FamilyState.Supervised)
					{
						int index = (fname == null)? 3: 4;
						int threads = Integer.parseInt(myargs[index]);
						boolean fastmode = Boolean.parseBoolean(myargs[index + 1]);
						boolean incore = Boolean.parseBoolean(myargs[index + 2]);
						ppicker = (fname == null)? new SupervisedParticlePicker(selfile, outputdir, threads, fastmode, incore): new SupervisedParticlePicker(selfile, outputdir, fname, threads, fastmode, incore);
						
					}

					else if (mode == FamilyState.Review)
					{
						String reviewfile = myargs[3];
						ppicker = (fname == null)? new ReviewParticlePicker(selfile, outputdir, reviewfile): new ReviewParticlePicker(selfile, outputdir, fname, reviewfile);

					}
					else if (mode == FamilyState.ReadOnly)
						ppicker = new ReadOnlyParticlePicker(selfile, outputdir);
					tp = new TrainingPickerJFrame(ppicker);
				}
				catch (Exception e)
				{
					JOptionPane.showMessageDialog(null, e.getMessage());
					
					
				}
			}
		});

	}

}