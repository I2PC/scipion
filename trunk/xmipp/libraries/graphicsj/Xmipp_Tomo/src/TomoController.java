/***************************************************************************
 *
 * @author: Jesus Cuenca (jcuenca@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

import ij.IJ;
import ij.WindowManager;
import ij.gui.GenericDialog;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.Method;
import javax.swing.SwingWorker;

// TODO: make TomoController extend AbstractController?

/**
 * - Why? 
 * I opted for Model-View-Controller as the overall UI paradigm, hence the Controller: here we
 * handle the interactive application workflow
 * @implements AdjustmentListener - slidebar events
 *   for buttons/commands, @see Command & TomoAction
 */
public class TomoController implements AdjustmentListener{
	private TomoWindow window;
	// modelToLoad: to avoid passing parameters in load BackgroundMethods
	private TomoData modelToLoad,model;
	
	// @see TomoController.playPause
	private boolean projectionsPlaying = false;
	// true if play follows a loop pattern (after last image, play first), false if pattern is ping-pong (after last, play last-1 and so on)
	private boolean playLoop=true;
	// play direction: forward (+1) or backward (-1)
	private int playDirection = 1;

	TomoController(TomoWindow window) {
		this.window = window;
	}

	/**
	 * - Why?
	 * For running Controller methods in the background
	 * @extends SwingWorker
	 */
	class BackgroundMethod extends SwingWorker<String, Object> {
		private Method backgroundMethod,doneMethod=null;
		private Object target;

		
		public BackgroundMethod(Method backgroundMethod, Object target) {
			this.backgroundMethod = backgroundMethod;
			this.target = target;
		}
		
		public BackgroundMethod(Object target,Method backgroundMethod, Method doneMethod) {
			this(backgroundMethod,target);
			this.doneMethod=doneMethod;
		}

		@Override
		public String doInBackground() {
			try {
				backgroundMethod.invoke(target);
			} catch (Exception ex) {
				Xmipp_Tomo.debug("BackgroundMethod.doInBackground()", ex);
			}
			return "";
		}

		@Override
		protected void done() {
			try{
				if(doneMethod!=null)
					doneMethod.invoke(target);
			} catch (Exception ex) {
				Xmipp_Tomo.debug("BackgroundMethod.done()", ex);
			}
		}
	} //end BackgroundMethod


	/**
	 * - Why?
	 *  hack to intercept ImageJ generic dialogs
	 */
	private class DialogCaptureThread extends Thread {
		private TomoWindow window;

		DialogCaptureThread(TomoWindow w) {
			window = w;
		}

		public void run() {
			try {
				sleep(TomoWindow.WAIT_FOR_DIALOG_SHOWS);
				window.captureDialog();
			} catch (Exception ex) {
				Xmipp_Tomo.debug(
						"DialogCaptureThread.run - unexpected exception", ex);
			}
		}
	} // end DialogCaptureThread

	
	public void playPause() {
		if (areProjectionsPlaying()) {
			// switch to pause
			setProjectionsPlaying(false);
			window.changeIcon(Command.PLAY.getId(), TomoWindow.PLAY_ICON);
		} else {
			try {
				new BackgroundMethod(getClass().getMethod("play"), this).execute();
			} catch (Exception ex) {
				Xmipp_Tomo.debug("TomoController.playPause() ", ex);
			}
		}
	}

	// this method must be public so Class.getMethod() finds it
	public void play() {
		setProjectionsPlaying(true);
		window.changeIcon(Command.PLAY.getId(), TomoWindow.PAUSE_ICON);
		while (areProjectionsPlaying()) {
			// 1..N
			int nextProjection = 0;

			if (isPlayLoop())
				nextProjection = getModel().getNextProjection();
			else {
				if ( getModel().getCurrentProjectionNumber() == getModel().getNumberOfProjections())
					setPlayDirection(-1);
				else if( getModel().getCurrentProjectionNumber() == 1)
					setPlayDirection(1);
				nextProjection = getModel().getCurrentProjectionNumber() + getPlayDirection();
			}
			getModel().setCurrentProjectionNumber(nextProjection);
			try {
				Thread.sleep(50);
			} catch (InterruptedException ex) {
			}
		}
	}

	public void changePlayMode() {
		setPlayLoop(!isPlayLoop());
	}
	
	public void printWorkflow() {
		Xmipp_Tomo.printWorkflow();
	}

	public void loadEM() {
		if (window.getLastCommandState() == Command.State.LOADING) {
			// if the button was pressed while loading, it means user wants to
			// cancel
			window.setLastCommandState(Command.State.CANCELED);
			window.getButton(Command.LOAD.getId()).setText(
					Command.LOAD.getLabel());
			Xmipp_Tomo.debug("loadEM - cancelled");
		} else {
			String path = FileDialog.openDialog("Load EM", window);
			if (path == null)
				return;
			loadImage(path);

		}
	}
	
	/**
	 * Remember to set modelToLoad before calling this method
	 */
	public void readImage() {
		String errorMessage = "";
		try {
			TiltSeriesIO.read(modelToLoad);
		} catch (FileNotFoundException ex) {
			Xmipp_Tomo.debug("ImportDataThread.run - ", ex);
		} catch (IOException ex) {
			Xmipp_Tomo.debug("ImportDataThread.run - Error opening file ",
					ex);
			errorMessage = "Error opening file";
		} catch (InterruptedException ex) {
			Xmipp_Tomo.debug(
					"ImportDataThread.run - Interrupted exception", ex);
		} catch (OutOfMemoryError err) {
			Xmipp_Tomo.debug("ImportDataThread.run - Out of memory"
					+ err.toString());
			errorMessage = "Out of memory";
		} catch (Exception ex) {
			Xmipp_Tomo.debug("ImportDataThread.run - unexpected exception",
					ex);
		} finally {
			if (modelToLoad.getNumberOfProjections() < 1) {
				modelToLoad.loadCanceled();
				TomoWindow.alert(errorMessage);
			}
		}
	}
	
	private void loadImage(String path){
		// free previous model (if any) ??
		setModel(null);
		modelToLoad=new TomoData(path);
		setModel(modelToLoad);
		window.getButton(Command.LOAD.getId()).setText("Cancel " + Command.LOAD.getLabel());
		window.setLastCommandState(Command.State.LOADING);

		try {
			new BackgroundMethod(this, getClass().getMethod("readImage"), getClass().getMethod("loadedEM")).execute();
			getModel().waitForFirstImage();
		} catch (Exception ex) {
			Xmipp_Tomo.debug("TomoController.loadImage() ", ex);
		}

		if (getModel().getNumberOfProjections() > 0) {
			window.setImagePlusWindow();

			window.setTitle(window.getTitle());
			window.addView();
			window.addControls();
			window.addUserAction(new UserAction(window.getWindowId(),
					"Load", getModel().getFileName()));
			window.enableButtonsAfterLoad();
		}
	}

	/**
	 * Post-actions (after all the projections are loaded)
	 */
	public void loadedEM() {
		if (window.getLastCommandState() == Command.State.LOADING) {
			window.setLastCommandState(Command.State.LOADED);
			window.getButton(Command.LOAD.getId()).setText(
					Command.LOAD.getLabel());
		}
	}

	public void alignManual() {

	}

	public void alignAuto() {

	}

	public void alignCorrelation() {

	}

	/**
	 * Apply this window's workflow 
	 */
	public void apply() {

		String path = FileDialog.saveDialog("Save...", window);
		if (path == null) {
			return;
		}
		
		// TODO: bug - modelToSave has wrong properties: path(not updated), width&height =  0, always resized
		TomoData originalModel = getModel(),modelToSave=originalModel;
		
		// Reopen the file only if the data was resized
		if (originalModel.isResized()) {
			modelToLoad = new TomoData(getModel().getFilePath());
			modelToSave=modelToLoad;
			// modelToLoad.addPropertyChangeListener(window);
			try {
				window.setLastCommandState(Command.State.RELOADING);
				new BackgroundMethod(this, getClass().getMethod("readImage"), null).execute();
				modelToLoad.waitForLastImage();
			} catch (Exception ex) {
				Xmipp_Tomo.debug("actionApply - unexpected exception", ex);
			}

			window.setLastCommandState(Command.State.LOADED);
		}
		
		// iterate through the user actions that make sense
		for (UserAction currentAction : Xmipp_Tomo.getWorkflow(window.getLastAction())) {
			if (currentAction.isNeededForFile()) {
				// Xmipp_Tomo.debug("Applying " + currentAction.toString());
				window.setStatus("Applying " + currentAction.getCommand());
				currentAction.getPlugin().run(modelToSave.getImage());
			}
		}

		// write
		saveFile(window, modelToSave, path);
	}

	public void gaussian() {
		Plugin plugin = new GaussianPlugin();
		window.setPlugin(plugin);
		runIjCmd(Command.GAUSSIAN, plugin);
	}

	public void measure() {
		Plugin plugin = new MeasurePlugin();
		window.setPlugin(plugin);
		runIjCmd(Command.MEASURE, plugin);
	}

	public void median() {
		Plugin plugin = new MedianPlugin();
		window.setPlugin(plugin);
		runIjCmd(Command.MEDIAN, plugin);
	}

	public void subBackground() {
		Plugin plugin = new BackgroundSubstractPlugin();
		window.setPlugin(plugin);
		runIjCmd(Command.SUB_BACKGROUND, plugin);
	}

	public void enhanceContrast() {
		Plugin plugin = new ContrastEnhancePlugin();
		window.setPlugin(plugin);
		runIjCmd(Command.ENHANCE_CONTRAST, plugin);
	}
	
	public void gammaCorrection() {
		Plugin plugin = new GammaCorrectionPlugin();
		window.setPlugin(plugin);
		runIjCmd(Command.GAMMA_CORRECTION, plugin);
	}
	
	public void bandpass() {
		Plugin plugin = new BandpassPlugin();
		window.setPlugin(plugin);
		runIjCmd(Command.BANDPASS, plugin);
	}
	
	public void histogramEqualization() {
		// convert to 8 bit
		getModel().convertTo(8);
		window.addUserAction(new UserAction(window.getWindowId(), "convertToByte"));
		
		// histogram equalization is an option of enhance contrast
		// @todo: it may be better to use a custom dialog, so user does not need to check "Histogram Equalization" and "Normalize all"
		Plugin plugin = new ContrastEnhancePlugin();
		window.setPlugin(plugin);
		runIjCmd(Command.ENHANCE_CONTRAST, plugin);
		
		// convert back to 32 bit
		window.protectWindow();
		getModel().convertTo(32);
		window.protectWindow();
		window.addUserAction(new UserAction(window.getWindowId(), "convertToFloat"));
		
		window.refreshImageCanvas();
	}
	
	// TODO: hotspotRemoval
	public void hotspotRemoval(){
		
	}
	
	// TODO: adjustbc
	public void adjustbc(){
		
	}
	
	public void crop(){
		window.protectWindow();
		
		runIjCmd("Crop");
		
		window.unprotectWindow();
		// TODO: adjust viewPanel to new image size
	}

	
	public void runIjCmd(String command) {

		window.setChangeSaved(false);

		WindowManager.setTempCurrentImage(getModel().getImage());
		try {
			IJ.run(command);
		} catch (Exception ex) {
			// catch java.lang.RuntimeException: Macro canceled
			Xmipp_Tomo.debug("actionRunIjCmd - action canceled");
			return;
		}
		window.refreshImageCanvas();

		window.setStatus("Done");

	}
	
	
	public void runIjCmd(Command command, Plugin plugin) {
		String label = command.getLabel();
		String cmd = plugin.getCommand();

		if (label != null)
			window.setStatus(label + " - started");

		window.setChangeSaved(false);

		if(plugin != null)
			(new Thread(new DialogCaptureThread(window))).start();

		runIjCmd(cmd);

		if(plugin != null)
			window.addUserAction(new UserAction(window.getWindowId(), cmd, window
				.getPlugin()));

		window.setPlugin(null);
	}

	// Allow user to just make a copy of a file
	public void save() {
		String path = FileDialog.saveDialog("Save...", window);
		if (path == null) {
			return;
		}
		// TODO: should include resize code (see apply() )
		saveFile(window, getModel(), path);
	}

	public void setTilt() {
		GenericDialog gd = new GenericDialog("Tilt angles");
		gd.addMessage("Please fill either start and end angles, or start and step");
		gd.addNumericField("Start angle", 0.0, 1);
		gd.addNumericField("End angle", 0.0, 1);
		gd.addNumericField("Step", 0.0, 1);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		double start = gd.getNextNumber();
		double end = gd.getNextNumber();
		double step = gd.getNextNumber();
		if (step == 0.0)
			getModel().setTiltAngles(start, end);
		else
			getModel().setTiltAnglesStep(start, step);
	}

	/**
	 * Ask user parameters for xmipp_xray_import and run it (by now, the program
	 * must be on $PATH)
	 */
	public void loadXray() {
		XrayImportDialog d = new XrayImportDialog("X-Ray import", window);
		// d.setup();
		d.showDialog();
		if (d.getStatus() == Xmipp_Tomo.ExitValues.OK) {
			String command = d.getCommand();
			// Xmipp_Tomo.debug("bash -c " + '"' + command + '"');
			Xmipp_Tomo.ExitValues error = Xmipp_Tomo.exec(command,false);
			if (error == Xmipp_Tomo.ExitValues.OK){
				String imagePath= d.getImagePath();
				loadImage(imagePath);
			}else
				Xmipp_Tomo.debug("Error (" + error + ") executing " + command);
		}

		
	}

	public void normalize(){
		getModel().normalize();
	}
	
	public void discardProjection(){
		if(getModel().isCurrentEnabled()){
			getModel().discardCurrentProjection();
			/* window.refreshImageCanvas();
			window.changeLabel(Command.DISCARD_PROJECTION.getId(), Command.UNDO_DISCARD_PROJECTION.getLabel()); */
		}else{
			getModel().enableCurrentProjection();
			/* window.refreshImageCanvas();
			window.changeLabel(Command.DISCARD_PROJECTION.getId(), Command.DISCARD_PROJECTION.getLabel()); */
		}
	}
	
	public void currentProjectionInfo(){
		Xmipp_Tomo.debug(getModel().getCurrentProjectionInfo());
	}
	
	/**
	 * Scrollbar events
	 * 
	 * @param e
	 */
	public synchronized void adjustmentValueChanged(AdjustmentEvent e) {
		// Xmipp_Tomo.debug("TomoWindow.adjustmentValueChanged"+e.getSource().toString());
		// TODO: AdjustmentListener - verify "projectionScrollbar")
		//if (e.getSource().toString().equals("projectionScrollbar")) {
			// scrollbar range is 1..N, like projections array
		if(getModel() != null)
			getModel().setCurrentProjectionNumber(e.getValue());
		// }
		notify();
	}

	public TomoData getModel() {
		return model;
	}

	public void setModel(TomoData model) {
		this.model = model;
		if(window != null && model != null)
			model.addPropertyChangeListener(window);
	}
	
	private boolean isPlayLoop() {
		return playLoop;
	}

	private void setPlayLoop(boolean playLoop) {
		this.playLoop = playLoop;
	}

	private int getPlayDirection() {
		return playDirection;
	}

	private void setPlayDirection(int playDirection) {
		this.playDirection = playDirection;
	}
	
	private boolean areProjectionsPlaying() {
		return projectionsPlaying;
	}

	private void setProjectionsPlaying(boolean playing) {
		this.projectionsPlaying = playing;
	}

	public void saveFile(TomoWindow tomoWindow, TomoData model, String path) {
		tomoWindow.setStatus("Saving...");
		model.setFile(path);
		TiltSeriesIO.write(model);
		tomoWindow.setStatus("Done");
		tomoWindow.setChangeSaved(true);
	}
	
	// TODO: xmipp_tomo_remove_fluctuations (preprocessing)
}
