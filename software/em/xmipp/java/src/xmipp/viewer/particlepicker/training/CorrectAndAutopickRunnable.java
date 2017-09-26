/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.particlepicker.training;

import java.util.logging.Level;
import xmipp.jni.MDRow;
import xmipp.jni.Particle;
import xmipp.jni.PickingClassifier;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.training.gui.SupervisedPickerJFrame;
import xmipp.viewer.particlepicker.training.model.MicrographState;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.viewer.particlepicker.training.model.SupervisedPickerMicrograph;

import javax.swing.*;

/**
 *
 * @author airen
 */
public class CorrectAndAutopickRunnable implements Runnable
	{

		private final MDRow[] manualRows;
		private final MDRow[] automaticRows;
        private Particle[] autopickRows;

		private SupervisedPickerMicrograph next;
		private SupervisedPickerJFrame frame;
        private final SupervisedParticlePicker picker;
        private final PickingClassifier classifier;
                

		public CorrectAndAutopickRunnable(SupervisedPickerJFrame frame, MDRow[] manualRows, MDRow[] automaticRows,
				SupervisedPickerMicrograph next)
		{
			this.frame = frame;
			this.manualRows = manualRows;
			this.automaticRows = automaticRows;
			this.next = next;
            this.picker = frame.getParticlePicker();
            this.classifier = (PickingClassifier)picker.getClassifier();
		}

		public void run()
		{
            try
            {
				classifier.correct(manualRows, automaticRows);
	
				if (picker.getMode() == Mode.Supervised && next.getState() == MicrographState.Supervised)
				{
					next.getAutomaticParticles().clear();
					next.setAutopickpercent(picker.getAutopickpercent());
					autopickRows = classifier.autopick(next.getFile(), next.getAutopickpercent());
					picker.loadParticles(next, autopickRows);
                    picker.saveData(next);
                    // Runs inside of the Swing UI thread
                    SwingUtilities.invokeLater(new Runnable() {
                        public void run() {
                            frame.setChanged(false);
                        }
                    });

				}


                // Runs inside of the Swing UI thread
                SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        // Select the worse particle as active
                        frame.getCanvas().refreshActive(null);
                        frame.getCanvas().repaint();
                        frame.getCanvas().setEnabled(true);
                        XmippWindowUtil.releaseGUI(frame.getRootPane());
                        frame.updateMicrographsModel(true);
                    }
                });

            }
            catch(Exception e)
            {

                ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
                SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        XmippDialog.showError(frame, "Classifier error");
                    }
                });

            }
		}
	}
