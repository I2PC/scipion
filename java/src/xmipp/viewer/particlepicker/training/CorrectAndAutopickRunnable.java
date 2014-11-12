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
import xmipp.viewer.particlepicker.training.gui.SupervisedParticlePickerJFrame;
import xmipp.viewer.particlepicker.training.model.MicrographState;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.viewer.particlepicker.training.model.SupervisedPickerMicrograph;

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
		private SupervisedParticlePickerJFrame frame;
                private final SupervisedParticlePicker picker;
                private final PickingClassifier classifier;
    private final SupervisedPickerMicrograph micrograph;
                

		public CorrectAndAutopickRunnable(SupervisedParticlePickerJFrame frame, MDRow[] manualRows, MDRow[] automaticRows,
				SupervisedPickerMicrograph next)
		{
			this.frame = frame;
			this.manualRows = manualRows;
			this.automaticRows = automaticRows;
			this.next = next;
                        this.micrograph = frame.getMicrograph();
                        this.picker = frame.getParticlePicker();
                        this.classifier = picker.getClassifier();
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
                                frame.setChanged(false);
			}
                        micrograph.resetParticlesRectangle();//after correct there is no need to keep this information for mic
			frame.getCanvas().repaint();
			frame.getCanvas().setEnabled(true);
			XmippWindowUtil.releaseGUI(frame.getRootPane());
                        frame.updateMicrographsModel(true);
                    }
                    catch(Exception e)
                    {
                        ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
                        XmippDialog.showError(frame, "Classifier error");
                    }
		}
	}
