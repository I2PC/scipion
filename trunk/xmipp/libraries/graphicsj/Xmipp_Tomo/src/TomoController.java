import ij.IJ;
import ij.WindowManager;
import ij.gui.GenericDialog;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.Method;

import javax.swing.SwingWorker;

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

/**
 * @author jcuenca Collection of GUI action methods
 */
public class TomoController {
	TomoWindow window;

	TomoController(TomoWindow window) {
		this.window = window;
	}

	class BackgroundMethod extends SwingWorker<String, Object> {
		Method method;
		Object target;

		public BackgroundMethod(Method m, Object target) {
			method = m;
			this.target = target;
		}

		@Override
		public String doInBackground() {
			try {
				method.invoke(target);
			} catch (Exception ex) {
				Xmipp_Tomo.debug("BackgroundMethod.doInBackground()");
			}
			return "";
		}

		@Override
		protected void done() {
		}
	}

	// ImportDataThread: load & display projections in parallel
	private class ImportDataThread extends Thread {
		private TomoData dataModel;
		private boolean resize = false;
		private TomoController controller;

		ImportDataThread(TomoData model, boolean resize,TomoController controller) {
			dataModel = model;
			this.resize = resize;
			this.controller = controller;
		}

		public void run() {
			try {
				new TiltSeriesIO(resize).read(dataModel);
				controller.loadedEM();
			} catch (FileNotFoundException ex) {
				Xmipp_Tomo.debug("ImportDataThread.run - " + ex.getMessage());
			} catch (IOException ex) {
				Xmipp_Tomo.debug("ImportDataThread.run - Error opening file");
			} catch (InterruptedException ex) {
				Xmipp_Tomo
						.debug("ImportDataThread.run - Interrupted exception");
			} catch (OutOfMemoryError err){
				Xmipp_Tomo.debug("ImportDataThread.run - Out of memory");
			} catch (Exception ex) {
				Xmipp_Tomo.debug("ImportDataThread.run - unexpected exception",
						ex);
			}
		}
	}

	// implement UI concurrency with private class that extends Thread
	// DialogCaptureThread: intercept generic dialogs
	private class DialogCaptureThread extends Thread {
		private TomoWindow window;

		DialogCaptureThread(TomoWindow w) {
			window = w;
		}

		public void run() {
			try {
				sleep(TomoWindow.WAIT_FOR_DIALOG);
				window.captureDialog();
			} catch (Exception ex) {
				Xmipp_Tomo.debug(
						"DialogCaptureThread.run - unexpected exception", ex);
			}
		}
	}

	// Play / pause automatic projection display in its own thread
	public void playPause() {
		if (window.isPlaying()) {
			// switch to pause
			window.setPlaying(false);
			window.changeIcon(Command.PLAY.getId(), TomoWindow.PLAY_ICON);
			// PLAY.getButton().setIcon(getIcon(PLAY_ICON));
		} else {
			try {
				new BackgroundMethod(getClass().getMethod("play"), this)
						.execute();
			} catch (Exception ex) {
				Xmipp_Tomo.debug("TomoController.playPause() "
						+ ex.getMessage());
			}
		}
	}

	// this method must be public so Class.getMethod() finds it
	public void play() {
		// switch to play
		window.setPlaying(true);
		window.changeIcon(Command.PLAY.getId(), TomoWindow.PAUSE_ICON);
		while (window.isPlaying()) {
			int nextProjection = (window.getModel().getCurrentProjection() + 1)
					% (window.getModel().getNumberOfProjections() + 1);
			window.getProjectionScrollbar().setValue(nextProjection);
			try {
				Thread.sleep(50);
			} catch (InterruptedException ex) {
			}
		}
	}

	public void printWorkflow() {
		Xmipp_Tomo.printWorkflow();
	}

	public void loadEM() {
		if (window.getLastCommandState() == Command.State.LOADING) {
			// if the button was pressed while loading, it means user wants to
			// cancel
			window.setLastCommandState(Command.State.CANCELED);
			window.getButton(Command.LOAD.getId()).setText(Command.LOAD.getLabel());
			Xmipp_Tomo.debug("loadEM - cancelled");
		} else {
			String path = XrayImportDialog.dialogOpen();
			if ((path == null) || ("".equals(path)))
				return;

			// free previous model (if any) ??
			window.setModel(null);
			window.setModel(new TomoData(path, window));
			window.getButton(Command.LOAD.getId()).setText("Cancel " + Command.LOAD.getLabel());
			window.setLastCommandState(Command.State.LOADING);
			
			try {
				// import data in one thread and hold this thread until the
				// first
				// projection is loaded
				
				(new Thread(new ImportDataThread(window.getModel(), true,this)))
						.start();
				window.getModel().waitForFirstImage();
			} catch (InterruptedException ex) {
				Xmipp_Tomo.debug("Xmipp_Tomo - Interrupted exception");
			}

			window.setImagePlusWindow();

			window.setTitle(window.getTitle());
			window.addView();
			window.addControls();
			window.addUserAction(new UserAction(window.getWindowId(), "Load",
					window.getModel().getFileName()));

			// enable buttons
			window.enableButtonsAfterLoad();

			// wait for last image in a different thread
			try {
				new BackgroundMethod(getClass().getMethod("loadedEM"), this)
						.execute();
			} catch (Exception ex) {
				Xmipp_Tomo.debug("TomoController.playPause() "
						+ ex.getMessage());
			}
		}
	}

	/**
	 * Post-actions (once all the images are loaded)
	 */
	public void loadedEM() {
		if(window.getLastCommandState() == Command.State.LOADING){
			window.setLastCommandState(Command.State.LOADED);
			window.getButton(Command.LOAD.getId()).setText(Command.LOAD.getLabel());
		}
	}

	public void alignManual() {

	}

	public void alignAuto() {

	}

	public void alignCorrelation() {

	}

	/**
	 * Apply this window's workflow to the file associated to this window
	 * (through the model)
	 */
	public void apply() {

		String path = window.dialogSave();
		if ("".equals(path)) {
			return;
		}

		// Reopen the file only if the data was resized
		TomoData originalModel = window.getModel();

		if (window.getModel().isResized()) {
			originalModel = new TomoData(window.getModel().getFilePath(),
					window);
			originalModel.addPropertyChangeListener(window);
			try {
				window.setLastCommandState(Command.State.RELOADING);
				(new Thread(new ImportDataThread(originalModel, false,this)))
						.start();
				originalModel.waitForLastImage();
			} catch (Exception ex) {
				Xmipp_Tomo.debug("actionApply - unexpected exception", ex);
			}

			window.setLastCommandState(Command.State.LOADED);
			// iterate through the user actions that make sense
			for (UserAction currentAction : Xmipp_Tomo.getWorkflow(window
					.getLastAction())) {
				if (currentAction.isNeededForFile()) {
					// Xmipp_Tomo.debug("Applying " + currentAction.toString());
					window.setStatus("Applying " + currentAction.getCommand());
					currentAction.getPlugin().run(originalModel.getImage());
				}
			}
		}

		// write
		window.saveFile(originalModel, path);
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

	public void runIjCmd(Command command, Plugin plugin) {
		String label = command.getLabel();
		String cmd = plugin.getCommand();

		if (label != null)
			window.setStatus(label + " - started");

		window.setChangeSaved(false);

		(new Thread(new DialogCaptureThread(window))).start();

		WindowManager.setTempCurrentImage(window.getModel().getImage());
		try {
			IJ.run(cmd);
		} catch (Exception ex) {
			// catch java.lang.RuntimeException: Macro canceled
			Xmipp_Tomo.debug("actionRunIjCmd - action canceled");
			return;
		}
		window.refreshImageCanvas();

		window.setStatus("Done");

		window.addUserAction(new UserAction(window.getWindowId(), cmd, window
				.getPlugin()));

		window.setPlugin(null);
	}

	public void save() {
		String path = window.dialogSave();
		if ((path == null) || ("".equals(path)))
			return;
		window.saveFile(window.getModel(), path);
	}

	public void setTilt() {
		GenericDialog gd = new GenericDialog("Tilt angles");
		gd
				.addMessage("Please fill either start and end angles, or start and step");
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
			// use start & end
			window.getModel().setTiltAngles(start, end);
		else
			window.getModel().setTiltAnglesStep(start, step);
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
			Xmipp_Tomo.ExitValues error = Xmipp_Tomo.exec(command);
			if (error != Xmipp_Tomo.ExitValues.OK)
				Xmipp_Tomo.debug("Error (" + error + ") executing " + command);
		}
	}
}
