package xmipp.particlepicker.tiltpair;



import javax.swing.SwingUtilities;

import xmipp.particlepicker.tiltpair.gui.TiltPairPickerJFrame;
import xmipp.particlepicker.tiltpair.model.TiltPairPicker;


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
					
				
				TiltPairPicker pppicker = new TiltPairPicker(selfile, outputdir);
				TiltPairPickerJFrame frame = new TiltPairPickerJFrame(pppicker);

			}
		});
	}

}