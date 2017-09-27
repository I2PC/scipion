/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.particlepicker.training;

import xmipp.jni.Classifier;
import xmipp.jni.Particle;
import xmipp.jni.PickingClassifier;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.training.gui.SupervisedPickerJFrame;
import xmipp.viewer.particlepicker.training.model.GenericClassifier;
import xmipp.viewer.particlepicker.training.model.MicrographState;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.viewer.particlepicker.training.model.SupervisedPickerMicrograph;

import javax.swing.*;

/**
 *
 * @author airen
 */
public class AutopickRunnable implements Runnable
	{

		private SupervisedPickerJFrame frame;
		private Particle[] autopickRows;
		private SupervisedPickerMicrograph micrograph;
        private final SupervisedParticlePicker picker;
        private final Classifier classifier;
                
		public AutopickRunnable(SupervisedPickerJFrame frame, SupervisedPickerMicrograph micrograph)
		{
			this.frame = frame;
			this.micrograph = micrograph;
            this.picker = frame.getParticlePicker();
            this.classifier = picker.getClassifier();
		}

		public void run()
		{
			micrograph.getAutomaticParticles().clear();
			if(classifier instanceof PickingClassifier)
			{
				micrograph.setState(MicrographState.Supervised);
				micrograph.setAutopickpercent(picker.getAutopickpercent());
				autopickRows = ((PickingClassifier)classifier).autopick(micrograph.getFile(), micrograph.getAutopickpercent());
				picker.loadParticles(micrograph, autopickRows);
	            picker.saveData(micrograph);
			}
			else 
			{
				((GenericClassifier)classifier).autopick(micrograph);
				picker.loadMicrographData(micrograph);
			}

            // Runs inside of the Swing UI thread
            SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    frame.setChanged(false);

                    // Select the worse particle as active
                    frame.getCanvas().refreshActive(null);

                    frame.getCanvas().repaint();
                    frame.getCanvas().setEnabled(true);
                    XmippWindowUtil.releaseGUI(frame.getRootPane());
                    frame.updateMicrographsModel();
                }
            });
		}

	}
