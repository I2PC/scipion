package xmipp.viewer.particlepicker.tiltpair;



import javax.swing.SwingUtilities;

import xmipp.viewer.particlepicker.tiltpair.gui.TiltPairPickerJFrame;
import xmipp.viewer.particlepicker.tiltpair.model.TiltPairPicker;
import xmipp.viewer.particlepicker.training.model.Mode;


class Main {
	// 0 --> input metadata
	// 1 --> output dir
	// 2 --> mode
	
	//On Supervised
	// 3 --> number of threads for supervised mode
	// 4 --> fast mode for supervised mode
	// 5 --> incore for supervised mode
	
	//On Review
	//3 -->particlesfile
	public static void main(String[] args) {

		final String[] myargs = args;
		
		
		SwingUtilities.invokeLater(new Runnable() {

			@Override
			public void run() {
				
				String selfile = myargs[0];
				String outputdir = myargs[1];
				Mode state = (myargs.length == 2)? Mode.Manual: Mode.getMode(myargs[2]);	
				
				TiltPairPicker pppicker = new TiltPairPicker(selfile, outputdir, state);
				TiltPairPickerJFrame frame = new TiltPairPickerJFrame(pppicker);

			}
		});
	}

}