package Xmipp;

import javax.swing.SwingUtilities;


import gui.ParticlePickerJFrame;
import ij.plugin.PlugIn;

public class Xmipp_Particle_Picker_ implements PlugIn {

	@Override
	public void run(String args){
			
			SwingUtilities.invokeLater(new Runnable() {
				
				@Override
				public void run() {
					try {
					//ParticlePickerJFrame frame = new ParticlePickerJFrame();
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
		});
	}

}
