package picker;



import javax.swing.SwingUtilities;

import picker.gui.ParticlePickerJFrame;
import picker.model.FamilyState;
import picker.model.ManualParticlePicker;
import picker.model.ParticlePicker;
import picker.model.ReviewParticlePicker;
import picker.model.SupervisedParticlePicker;


class Main {
	// 0 --> input metadata
	// 1 --> output dir
	// 2 --> mode
	
	//On Supervised
	// 3 --> number of threads for supervised mode
	// 4 --> fast mode for supervised mode
	// 5 --> incore for supervised mode
	
	//On Review
	//3 -->external auto dir
	public static void main(String[] args) {

		final String[] myargs = args;
		
		
		SwingUtilities.invokeLater(new Runnable() {

			@Override
			public void run() {
				ParticlePicker ppicker = null;
				String selfile = myargs[0];
				String outputdir = myargs[1];
				FamilyState mode = FamilyState.getFamilyState(myargs[2]);
		
				if(mode == FamilyState.Manual)
					ppicker = new ManualParticlePicker(selfile, outputdir, mode);
				
				if (mode == FamilyState.Supervised) {
					int threads = Integer.parseInt(myargs[3]);
					boolean fastmode = Boolean.parseBoolean(myargs[4]);
					boolean incore = Boolean.parseBoolean(myargs[5]);
					ppicker = new SupervisedParticlePicker(selfile, outputdir, threads, fastmode, incore);
				}
				
				
				else if(mode == FamilyState.Review)
				{
					String reviewfile = myargs[3];
					ppicker = new ReviewParticlePicker(selfile, outputdir, reviewfile);
				}
				new ParticlePickerJFrame(ppicker);

			}
		});
	}

}