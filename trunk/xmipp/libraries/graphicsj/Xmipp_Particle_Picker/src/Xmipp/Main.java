package Xmipp;

import gui.ParticlePickerJFrame;

import javax.swing.SwingUtilities;

import model.ParticlePicker;

class Main {
	// 0 --> input metadata
	// 1 --> output dir
	// 2 --> number of threads for automatic picking
	// 3 --> fast mode for automatic picking
	// 4 --> incore for automatic picking
	public static void main(String[] args) {
		String mgselfile = args[0];
		String outputdir = args[1];
		ParticlePicker.setMicrographsSelFile(mgselfile);
		ParticlePicker.setOutputDir(outputdir);

		if (args.length == 2) {
			ParticlePicker.setIsAuto(false);
			SwingUtilities.invokeLater(new Runnable() {

				@Override
				public void run() {
					new ParticlePickerJFrame();

				}
			});
		} else {
			int threads = Integer.parseInt(args[2]);
			boolean fastMode = Boolean.parseBoolean(args[3]);
			boolean incore = Boolean.parseBoolean(args[4]);
			ParticlePicker.setThreads(threads);
			ParticlePicker.setFastMode(fastMode);
			ParticlePicker.setIncore(incore);
			ParticlePicker.setIsAuto(true);
			SwingUtilities.invokeLater(new Runnable() {

				@Override
				public void run() {
					new ParticlePickerJFrame();

				}
			});
		}

	}

}