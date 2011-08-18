package Xmipp;

import javax.swing.SwingUtilities;

import ppicker.XmippParticlePickerJFrame;

import ij.plugin.PlugIn;

public class Xmipp_Particle_Picker_ implements PlugIn {

	@Override
	public void run(String args){
			
			SwingUtilities.invokeLater(new Runnable() {
				
				@Override
				public void run() {
					try {
					XmippParticlePickerJFrame frame = new XmippParticlePickerJFrame();
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
		});
	}

}
