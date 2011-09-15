package tiltpairpicker;



import javax.swing.SwingUtilities;


import tiltpairpicker.gui.ParticlePairPickerJFrame;
import tiltpairpicker.model.ParticlePairPicker;
import trainingpicker.gui.ParticlePickerJFrame;
import trainingpicker.model.FamilyState;
import trainingpicker.model.ManualParticlePicker;
import trainingpicker.model.ParticlePicker;
import trainingpicker.model.ReviewParticlePicker;
import trainingpicker.model.SupervisedParticlePicker;


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
				
				String selfile = myargs[0];
				String outputdir = myargs[1];
				
				ParticlePairPicker pppicker = new ParticlePairPicker(selfile, outputdir);
				new ParticlePairPickerJFrame(pppicker);

			}
		});
	}

}