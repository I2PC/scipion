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
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.process.ImageProcessor;

import java.awt.Dimension;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.Method;
import java.util.Enumeration;

import javax.swing.SwingWorker;
import javax.swing.tree.DefaultMutableTreeNode;

import xmipp.CastWriteMode;
import xmipp.ImageDouble;
import xmipp.ImageWriteMode;
import xmipp.MDLabel;
import xmipp.MetaData;

// TODO: make TomoController extend AbstractController?
/**
 * 	// TODO: write operations (below)
 *  the "problem": most files are single files (for instance, an mrc). but selfiles imply 2 files:
 *  the selfile and the stack file. So, what file naming policy to use?
 // Sel: One option would be to ask for sel and stack file paths (defaults to same paths, 
 //      ask if user wants to reuse the image file, or overwrite it, or use a different name...)
 //      right now we reuse the image filename
 //      Other option is use same name and offer only to choose the stack type (stk or mrc)
 * Let's start with second one (choose type)
 // Stack: ask only for stack file path (same overwrite warning), then save both stack and tlt [OK]
 */

/**
 * - Why? I opted for Model-View-Controller as the overall UI paradigm, hence
 * the Controller: here we handle the interactive application workflow
 * 
 * @implements AdjustmentListener - slidebar events for buttons/commands, @see
 *             Command & TomoAction
 */
public class TomoController implements AdjustmentListener {
	private TomoWindow window;
	// modelToLoad: to avoid passing parameters in load BackgroundMethods
	private TomoData modelToLoad, model, modelToSave;
	private StackModel stackModel;

	public StackModel getStackModel() {
		return stackModel;
	}

	public void setStackModel(StackModel stackModel) {
		this.stackModel = stackModel;
	}

	// Workflow model
	private Workflow workflow;

	public Workflow getWorkflow() {
		return workflow;
	}

	public TomoData getModelToSave() {
		return modelToSave;
	}

	public void setModelToSave(TomoData modelToSave) {
		this.modelToSave = modelToSave;
	}

	private String destinationPath;

	public String getDestinationPath() {
		return destinationPath;
	}

	public void setDestinationPath(String destinationPath) {
		this.destinationPath = destinationPath;
	}

	// @see TomoController.playPause
	private boolean projectionsPlaying = false;
	// true if play follows a loop pattern (after last image, play first), false
	// if pattern is ping-pong (after last, play last-1 and so on)
	private boolean playLoop = true;
	// play direction: forward (+1) or backward (-1)
	private int playDirection = 1;

	// if an image is bigger than this threshold, resize it to this size
	// This optimization makes sense since for interactive use no more detail is
	// needed
	public static Dimension resizeThreshold = new Dimension(400, 400);

	TomoController(TomoWindow window) {
		this.window = window;
		workflow = new Workflow();
	}

	/**
	 * - Why? For running Controller methods in the background
	 * 
	 * @extends SwingWorker
	 */
	// TODO: allow parameters to backgroundmethods
	class BackgroundMethod extends SwingWorker<String, Object> {
		private Method backgroundMethod, doneMethod = null;
		private Object target;

		public BackgroundMethod(Method backgroundMethod, Object target) {
			this.backgroundMethod = backgroundMethod;
			this.target = target;
		}

		public BackgroundMethod(Object target, Method backgroundMethod,
				Method doneMethod) {
			this(backgroundMethod, target);
			this.doneMethod = doneMethod;
		}

		@Override
		public String doInBackground() {
			try {
				backgroundMethod.invoke(target);
			} catch (Exception ex) {
				Logger.debug("BackgroundMethod.doInBackground()", ex);
			}
			return "";
		}

		@Override
		protected void done() {
			try {
				if (doneMethod != null)
					doneMethod.invoke(target);
			} catch (Exception ex) {
				Logger.debug("BackgroundMethod.done()", ex);
			}
		}
	} // end BackgroundMethod

	/**
	 * - Why? hack to intercept ImageJ generic dialogs
	 */
	// TODO: pass the plugin as a parameter (instead of a global variable)
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
				Logger.debug("DialogCaptureThread.run - unexpected exception",
						ex);
			}
		}
	} // end DialogCaptureThread

	public void playPause() {
		if (areProjectionsPlaying()) {
			// switch to pause
			setProjectionsPlaying(false);
			window.changeIcon(XmippTomoCommands.PLAY.getId(),
					TomoWindow.PLAY_ICON);
		} else {
			try {
				new BackgroundMethod(getClass().getMethod("play"), this)
						.execute();
			} catch (Exception ex) {
				Logger.debug("TomoController.playPause() ", ex);
			}
		}
	}

	// this method must be public so Class.getMethod() finds it
	public void play() {
		setProjectionsPlaying(true);
		window
				.changeIcon(XmippTomoCommands.PLAY.getId(),
						TomoWindow.PAUSE_ICON);
		while (areProjectionsPlaying()) {
			// 1..N
			int nextProjection = 0;

			if (isPlayLoop())
				nextProjection = getStackModel().getNextProjection();
			else {
				if (getStackModel().getCurrentProjectionNumber() == getStackModel()
						.getNumberOfProjections())
					setPlayDirection(-1);
				else if (getStackModel().getCurrentProjectionNumber() == 1)
					setPlayDirection(1);
				nextProjection = getStackModel().getCurrentProjectionNumber()
						+ getPlayDirection();
			}
			getStackModel().setCurrentProjectionNumber(nextProjection);
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
		Logger.debug(getWorkflow().toString());
	}

	// TODO: reset workflow after loading a stack
	public void loadEM() {
		if (window.getLastCommandState() == Command.State.LOADING) {
			// if the button was pressed while loading, it means user wants to
			// cancel
			window.setLastCommandState(Command.State.CANCELED);
			window.getButton(XmippTomoCommands.LOAD.getId()).setText(
					XmippTomoCommands.LOAD.getLabel());
			Logger.debug("loadEM - cancelled");
		} else {
			String path = TomoFileDialog.openDialog("Load EM", window);
			if (path == null)
				return;
			// free previous model (if any) ??
			// setModel(null);
			// TODO: old tomodata model will be unnecesary with StackModel, hence remove next line when that moment arrives
			setModel(new TomoData(path));
			setStackModel(new StackModel(getWorkflow()));
			getStackModel().addPropertyChangeListener(window);
			/*
			 * setModelToLoad(new TomoData(path)); try { new
			 * BackgroundMethod(this, getClass().getMethod("loadEMBackground"),
			 * getClass().getMethod("loadedEM")).execute();
			 * getModelToLoad().waitForFirstImage(); } catch (Exception ex) {
			 * Logger.debug("TomoController.loadEM() ", ex); }
			 * 
			 * setModel(getModelToLoad());
			 * 
			 * if (getStackModel().getNumberOfProjections() > 0) {
			 * window.setImagePlusWindow();
			 * 
			 * window.setTitle(window.getTitle()); window.addView();
			 * window.addControls(); window.addUserAction(new
			 * UserAction(window.getWindowId(), "Load", "Load",
			 * getStackModel().getFileName())); window.enableButtonsAfterLoad();
			 * }
			 */

			UserAction ua = new UserAction(window.getWindowId(), "Load",
					"Load", path);
			ua.getIoDetails().setInputFilePath(path);
			// TODO: move addUserAction to workflow model (instead of the
			// window). Right now it's inside window because
			// insertion depends on the workflowview selected node
			window.addUserAction(ua);
			getWorkflow().setSelectedUserAction(ua);
			window.addView();
			window.addControls();
			window.enableButtonsAfterLoad();

		}
	}

	/**
	 * TODO: new approach: load only single projections on demand and store them
	 * in a cache (with prefetch) - see Xmipp_Tomo.java
	 */
	/**
	 *  @deprecated
	 */
	public void loadEMBackground() {
		// setModel(getModelToLoad());
		window.getButton(XmippTomoCommands.LOAD.getId()).setText(
				"Cancel " + XmippTomoCommands.LOAD.getLabel());
		window.setLastCommandState(Command.State.LOADING);

		String errorMessage = "";
		try {
			String absolutePath = getModelToLoad().getFilePath();

			// return if user canceled open dialog
			if ((absolutePath == null) || (absolutePath.equals("")))
				return;

			getModelToLoad().readMetadata(absolutePath);

			long ids[] = getModelToLoad().getStackIds();
			for (long projectionId : ids) {
				String filePath = getModelToLoad().getFilePath(projectionId);
				if (filePath != null) {
					ImageDouble image = new ImageDouble();
					image.readSlice(filePath);
					ImagePlus imagePlus = Converter.convertToImagePlus(image);
					getModelToLoad().setOriginalHeight(image.getYsize());
					getModelToLoad().setOriginalWidth(image.getXsize());

					if (shouldResize(image.getXsize(), image.getYsize())) {
						// resize/scale - use an aux image processor for all
						// projections
						ImageProcessor ip = imagePlus.getProcessor();
						ImageProcessor ipResized = ip.resize(
								resizeThreshold.width, resizeThreshold.height);
						imagePlus = new ImagePlus(filePath, ipResized);
						getModelToLoad().addProjection(imagePlus);
						getModelToLoad().setResized(true);
					} else
						getModelToLoad().addProjection(imagePlus);
				}
			}
			getModelToLoad().lastImageLoaded();

		} catch (FileNotFoundException ex) {
			Logger.debug("loadEMBackground - ", ex);
		} catch (IOException ex) {
			Logger.debug("loadEMBackground - Error opening file ", ex);
			errorMessage = "Error opening file";
		} catch (OutOfMemoryError err) {
			Logger.debug("loadEMBackground - Out of memory" + err.toString());
			errorMessage = "Out of memory";
		} catch (Exception ex) {
			Logger.debug("loadEMBackground - unexpected exception", ex);
		} finally {
			if (getModelToLoad().getNumberOfProjections() < 1) {
				getModelToLoad().loadCanceled();
				TomoWindow.alert(errorMessage);
			}
		}
	}

	// TODO: shouldResize - maybe there's a better place for this method?
	public static boolean shouldResize(int width, int height) {

		return (width > resizeThreshold.width)
				|| (height > resizeThreshold.height);
	}

	/**
	 * Post-actions (after all the projections are loaded)
	 */
	/**
	 * @deprecated
	 */
	public void loadedEM() {
		if (window.getLastCommandState() == Command.State.LOADING) {
			window.setLastCommandState(Command.State.LOADED);
			window.getButton(XmippTomoCommands.LOAD.getId()).setText(
					XmippTomoCommands.LOAD.getLabel());
		}
	}

	public void alignManual() {

	}

	public void alignAuto() {

	}

	public void alignCorrelation() {

	}

	/**
	 * Apply this window's workflow TODO: apply the selected "workflow path" in
	 * the workflow view
	 */
	/**
	 * @deprecated
	 */
	// TODO: in write actions, handle writeSel
	// TODO: in write actions, if the file exists delete it first (to avoid
	// overwrite problems)
	// TODO: progress bar while applying filters - or at least update status bar
	// TODO: update file name (model) in window title (once saved)
	// TODO: enhance contrast seems to not be applied...
	public void apply() {

		String destinationPath = TomoFileDialog.saveDialog("Save stack...",
				window);
		if (destinationPath == null) {
			return;
		}
		setDestinationPath(destinationPath);
		try {
			new BackgroundMethod(this, getClass().getMethod("applyBackground"),
					null).execute();
		} catch (Exception ex) {
			Logger.debug("TomoController.apply() ", ex);
		}

	}

/**
 * @deprecated
 * @param image
 * @return
 */
	private ImagePlus applyWorkflowTo(ImagePlus image) {
		// iterate through the user actions that make sense
		Enumeration e = getWorkflow().getRoot().breadthFirstEnumeration();
		while (e.hasMoreElements()) {
			UserAction currentAction = (UserAction) (((DefaultMutableTreeNode) e
					.nextElement()).getUserObject());
			if (currentAction.isNeededForFile()) {
				// Xmipp_Tomo.debug("Applying " + currentAction.toString());
				// window.setStatus("Applying " + currentAction.getCommand());
				currentAction.getPlugin().run(image);
			}
		}
		return image;
	}

	/**
	 * @deprecated
	 */
	public void applyBackground() {
		// TODO: bug - modelToSave has wrong properties: path(not updated),
		// width&height = 0, always resized
		String sourcePath = getModel().getFilePath();
		try {
			if ((sourcePath == null) || (sourcePath.equals("")))
				throw new IOException("Empty path");

			TomoData modelToSave = getModel();
			modelToSave.readMetadata(sourcePath);

			long ids[] = modelToSave.getStackIds();
			for (long projectionId : ids) {
				window.setStatus("Applying workflow and saving..."
						+ projectionId + "/"
						+ getModel().getNumberOfProjections());
				String sourceProjectionPath = modelToSave
						.getFilePath(projectionId);
				if (sourceProjectionPath != null) {
					ImageDouble image = new ImageDouble();
					image.readSlice(sourceProjectionPath);
					ImagePlus ip = Converter.convertToImagePlus(image);
					ip = applyWorkflowTo(ip);
					image = Converter.convertToImageDouble(ip);
					// the recommended way is using the @ notation,
					// TODO: get the @ path with the help of the Filename xmipp
					// class
					String sliceDestPath = String.valueOf(projectionId) + "@"
							+ getDestinationPath();
					// Xmipp_Tomo.debug(sliceDestPath);
					image.write(sliceDestPath);
					// image.write(destinationPath, (int)projectionId, false,
					// ImageWriteMode.WRITE_APPEND, CastWriteMode.CW_CAST);
				}
			}
		} catch (IOException ex) {
			Logger.debug("apply - Error opening file ", ex);
		} catch (Exception ex) {
			Logger.debug("apply ", ex);
		}
		window.setStatus("Done");
	}

	private TomoData reOpenIfResized(TomoData originalModel) {
		TomoData modelToSave = originalModel;
		if (originalModel.isResized()) {
			modelToLoad = new TomoData(getModel().getFilePath());
			modelToSave = modelToLoad;
			// modelToLoad.addPropertyChangeListener(window);
			try {
				window.setLastCommandState(Command.State.RELOADING);
				new BackgroundMethod(this, getClass().getMethod("readImage"),
						null).execute();
				modelToLoad.waitForLastImage();
			} catch (Exception ex) {
				Logger.debug("actionApply - unexpected exception", ex);
			}

			window.setLastCommandState(Command.State.LOADED);
		}
		return modelToSave;
	}

	public void gaussian() {
		Plugin plugin = new GaussianPlugin();
		window.setPlugin(plugin);
		runIjCmd(XmippTomoCommands.GAUSSIAN, plugin);
	}

	public void measure() {
		Plugin plugin = new MeasurePlugin();
		window.setPlugin(plugin);
		runIjCmd(XmippTomoCommands.MEASURE, plugin);
	}

	public void median() {
		Plugin plugin = new MedianPlugin();
		window.setPlugin(plugin);
		runIjCmd(XmippTomoCommands.MEDIAN, plugin);
	}

	public void subBackground() {
		Plugin plugin = new BackgroundSubstractPlugin();
		window.setPlugin(plugin);
		runIjCmd(XmippTomoCommands.SUB_BACKGROUND, plugin);
	}

	public void enhanceContrast() {
		Plugin plugin = new ContrastEnhancePlugin();
		window.setPlugin(plugin);
		runIjCmd(XmippTomoCommands.ENHANCE_CONTRAST, plugin);
	}

	public void gammaCorrection() {
		Plugin plugin = new GammaCorrectionPlugin();
		window.setPlugin(plugin);
		runIjCmd(XmippTomoCommands.GAMMA_CORRECTION, plugin);
	}

	public void bandpass() {
		Plugin plugin = new BandpassPlugin();
		window.setPlugin(plugin);
		runIjCmd(XmippTomoCommands.BANDPASS, plugin);
	}

	// TODO: image viewer pops out
	public void histogramEqualization() {
		// convert to 8 bit
		getModel().convertTo(8);
		window.addUserAction(new UserAction(window.getWindowId(),
				"convertToByte"));

		// histogram equalization is an option of enhance contrast
		// TODO: it may be better to use a custom dialog, so user does not need
		// to check "Histogram Equalization" and "Normalize all"
		Plugin plugin = new ContrastEnhancePlugin();
		window.setPlugin(plugin);
		runIjCmd(XmippTomoCommands.ENHANCE_CONTRAST, plugin);

		// convert back to 32 bit
		window.protectWindow();
		getModel().convertTo(32);
		window.protectWindow();
		window.addUserAction(new UserAction(window.getWindowId(),
				"convertToFloat"));

		window.refreshImageCanvas();
	}

	// TODO: hotspotRemoval
	public void hotspotRemoval() {

	}

	// TODO: adjustbc
	public void adjustbc() {

	}

	public void crop() {
		window.protectWindow();

		runIjCmd("Crop");

		window.unprotectWindow();
		// TODO: adjust viewPanel to new image size
	}

	public void runIjCmd(String command) {

		window.setChangeSaved(false);

		ImagePlus originalImage = getStackModel().getCurrentImage().duplicate();
		
		WindowManager.setTempCurrentImage(getStackModel().getCurrentImage());
		try {
			IJ.run(command);
		} catch (Exception ex) {
			// catch java.lang.RuntimeException: Macro canceled
			Logger.debug("actionRunIjCmd - action canceled");
			return;
		}
		// the preview modifies the image in the cache...
		getStackModel().updateCurrentImage(originalImage);
		window.refreshImageCanvas();

		// window.setStatus("Done");

	}

	public void runIjCmd(Command command, Plugin plugin) {
		String label = command.getLabel();
		String cmd = plugin.getCommand();

		if (label != null)
			window.setStatus(label + " - started");

		window.setChangeSaved(false);

		if (plugin != null)
			(new Thread(new DialogCaptureThread(window))).start();

		runIjCmd(cmd);

		if (plugin != null) {
			UserAction currentAction = getWorkflow().getSelectedUserAction();
			UserActionIO currentIo = currentAction.getIoDetails();
			
			UserAction newAction = new UserAction(window.getWindowId(), label,cmd,plugin);
			window.addUserAction(newAction);

			// TODO: find a better place for this  2 lines (maybe UserActionIO?)
			String newActionDir = getWorkflow().getWorkingDir() + "/" + newAction.getWorkingDir();
			new File(newActionDir).mkdirs();

			String output = newActionDir + "/" + currentIo.getInputFileName();
			Logger.debug("output file: " + output);

			UserActionIO uaio = new UserActionIO();
			uaio.setInputFilePath(output);
			newAction.setIoDetails(uaio);

			// TODO: run applyIJPlugin in the background
			/*
			 * try { new BackgroundMethod(this,
			 * getClass().getMethod("applyIJPlugin"), null).execute(); } catch
			 * (Exception ex) { Logger.debug("runIjCmd() ", ex); }
			 */
			applyIJPlugin(cmd, currentAction, newAction);

		}

		window.setPlugin(null);
	}

	// TODO:-current- applyIJplugin & runIjcmd
	private void applyIJPlugin(String command, UserAction currentAction,UserAction newAction) {
		String outputFile = newAction.getIoDetails().getInputFilePath();
		Plugin plugin = newAction.getPlugin();
		int numberOfProjections = currentAction.getIoDetails().getNumberOfProjections();
		for (int i=1; i < numberOfProjections; i++) {
			long projectionId = currentAction.getIoDetails().getProjectionId(i);
			newAction.setProgress("Converting..." + projectionId + "/"	+ numberOfProjections);
			String sourceProjectionPath = currentAction.getIoDetails().getFilePath(i);
			try {
				if (sourceProjectionPath != null) {
					ImageDouble image = new ImageDouble();
					image.readSlice(sourceProjectionPath);

					ImagePlus ip = Converter.convertToImagePlus(image);
					plugin.run(ip);
					
					image = Converter.convertToImageDouble(ip);
					
					String outputFullPath = String.valueOf(projectionId) + "@"+ outputFile;
					Logger.debug(outputFullPath);
					image.write(outputFullPath);
					
					getStackModel().updateImage(outputFullPath, image, ip);
				}
			} catch (Exception ex) {
				Logger.debug("Problem applying plugin to img " + projectionId,ex);
			}
		}
	}

	public void convert() {
		String destinationPath = TomoFileDialog
				.saveDialog("Convert...", window);
		if (destinationPath == null) {
			return;
		}
		setDestinationPath(destinationPath);
		/*
		 * TomoData originalModel = getModel(); TomoData
		 * modelToSave=reOpenIfResized(originalModel); // write saveFile(window,
		 * modelToSave, path);
		 */

		setModelToSave(getModel());

		try {
			new BackgroundMethod(this, getClass()
					.getMethod("convertBackground"), null).execute();
		} catch (Exception ex) {
			Logger.debug("TomoController.convert() ", ex);
		}
	}

	public void convertBackground() {
		try {
			if (getModel().isResized()) {
				String sourcePath = getModel().getFilePath();

				if ((sourcePath == null) || (sourcePath.equals("")))
					throw new IOException("Empty path");

				getModelToSave().readMetadata(sourcePath);

				long ids[] = getModelToSave().getStackIds();
				for (long projectionId : ids) {
					window.setStatus("Converting..." + projectionId + "/"
							+ getModel().getNumberOfProjections());
					String sourceProjectionPath = getModelToSave().getFilePath(
							projectionId);
					if (sourceProjectionPath != null) {
						ImageDouble image = new ImageDouble();
						image.readSlice(sourceProjectionPath);
						// the recommended way is using the @ notation,
						// TODO: get the @ path with the help of the Filename
						// xmipp class
						String sliceDestPath = String.valueOf(projectionId)
								+ "@" + destinationPath;
						Logger.debug(sliceDestPath);
						image.write(sliceDestPath);
						// image.write(destinationPath, (int)projectionId,
						// false, ImageWriteMode.WRITE_APPEND,
						// CastWriteMode.CW_CAST);
					}
				}
			} else {
				ImageDouble img = Converter.convertToImageDouble(getModel()
						.getImage());
				img.setFilename(getDestinationPath());
				// @see rwSpider.cpp -- WRITE_OVERWRITE required for maxim to be
				// updated (readSpider gets nDim from maxim)
				img.write(getDestinationPath(), 0, true,
						ImageWriteMode.WRITE_OVERWRITE, CastWriteMode.CW_CAST);
			}
		} catch (IOException ex) {
			Logger.debug("convert - Error opening file ", ex);
		} catch (Exception ex) {
			Logger.debug("convert ", ex);
		}

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
			getModel().setTiltAngles(start, end);
		else
			getModel().setTiltAnglesStep(start, step);
	}

	/**
	 * Ask user parameters for xmipp_xray_import and run it (by now, the program
	 * must be on $PATH)
	 */
	// TODO: progress bar while running xmipp_xray_import
	public void loadXray() {
		XrayImportDialog d = new XrayImportDialog("X-Ray import", window);
		// d.setup();
		d.showDialog();
		if (d.getStatus().equals(ExitValue.OK)) {
			String command = d.getCommand();
			// Xmipp_Tomo.debug("bash -c " + '"' + command + '"');
			ExitValue error = Xmipp_Tomo.exec(command, false);
			if (error == ExitValue.OK) {
				String imagePath = d.getImagePath();
				// loadImage(imagePath);
				setModelToLoad(new TomoData(imagePath));
				try {
					new BackgroundMethod(this, getClass().getMethod(
							"loadEMBackground"), getClass().getMethod(
							"loadedEM")).execute();
					getModel().waitForFirstImage();
				} catch (Exception ex) {
					Logger.debug("TomoController.loadEM() ", ex);
				}
			} else
				Logger.debug("Error (" + error + ") executing " + command);
		}

	}

	public void normalize() {
		getModel().normalize();
	}

	public void discardProjection() {
		if (getModel().isCurrentEnabled()) {
			getModel().discardCurrentProjection();
			/*
			 * window.refreshImageCanvas();
			 * window.changeLabel(Command.DISCARD_PROJECTION.getId(),
			 * Command.UNDO_DISCARD_PROJECTION.getLabel());
			 */
		} else {
			getModel().enableCurrentProjection();
			/*
			 * window.refreshImageCanvas();
			 * window.changeLabel(Command.DISCARD_PROJECTION.getId(),
			 * Command.DISCARD_PROJECTION.getLabel());
			 */
		}
	}

	public void currentProjectionInfo() {
		Logger.debug(getModel().getCurrentProjectionInfo());
	}

	/**
	 * Scrollbar events
	 * 
	 * @param e
	 */
	public synchronized void adjustmentValueChanged(AdjustmentEvent e) {
		// Xmipp_Tomo.debug("TomoWindow.adjustmentValueChanged"+e.getSource().toString());
		// TODO: AdjustmentListener - verify "projectionScrollbar")
		// if (e.getSource().toString().equals("projectionScrollbar")) {
		// scrollbar range is 1..N, like projections array
		if (getStackModel() != null)
			getStackModel().setCurrentProjectionNumber(e.getValue());
		// }
		notify();
	}

	public TomoData getModel() {
		return model;
	}

	public void setModel(TomoData model) {
		this.model = model;
		/*
		 * if(window != null && model != null)
		 * model.addPropertyChangeListener(window);
		 */
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

	public TomoData getModelToLoad() {
		return modelToLoad;
	}

	public void setModelToLoad(TomoData modelToLoad) {
		this.modelToLoad = modelToLoad;
	}

	// TODO: xmipp_tomo_remove_fluctuations (preprocessing)
}
