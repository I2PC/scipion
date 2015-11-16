/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.particlepicker.training;

import java.awt.Rectangle;
import java.util.logging.Level;
import xmipp.jni.MDRow;
import xmipp.jni.Particle;
import xmipp.jni.PickingClassifier;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.training.gui.SupervisedPickerJFrame;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.viewer.particlepicker.training.model.SupervisedPickerMicrograph;

/**
 *
 * @author airen
 */
public class TrainRunnable implements Runnable
	{

		protected SupervisedPickerJFrame frame;
		protected MDRow[] trainInput;
		protected Rectangle rectangle;
		protected Particle[] autopickRows;
		protected SupervisedPickerMicrograph micrograph;
		protected final SupervisedParticlePicker picker;
		protected final PickingClassifier classifier;

		public TrainRunnable(SupervisedPickerJFrame frame, MDRow[] trainInput, SupervisedPickerMicrograph trainmic)
		{
			this.frame = frame;
			this.trainInput = trainInput;
            this.micrograph = frame.getMicrograph();
            rectangle = trainmic.getRectangle();
            if(trainmic != micrograph)
            	trainmic.resetParticlesRectangle();
            this.picker = frame.getParticlePicker();
            this.classifier = (PickingClassifier)picker.getClassifier();
		}

		public void run()
		{
			try
			{
				classifier.train(trainInput, (int) rectangle.getX(), (int) rectangle.getY(), (int) rectangle.getWidth(), (int) rectangle.getHeight());// should remove training
				micrograph.setAutopickpercent(picker.getAutopickpercent());
				autopickRows = classifier.autopick(micrograph.getFile(), micrograph.getAutopickpercent());
				picker.loadParticles(micrograph, autopickRows);
                picker.saveData(micrograph);
				XmippWindowUtil.releaseGUI(frame.getRootPane());
				frame.getCanvas().setEnabled(true);
				frame.getCanvas().repaint();
				frame.updateMicrographsModel();
			}
			catch (Exception e)
			{
				ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
                                XmippDialog.showError(frame, "Classifier error");
			}
		}

		
	}