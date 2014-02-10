/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package xmipp.viewer.windows;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageWindow;
import ij.process.StackConverter;
import ij3d.Content;
import ij3d.Image3DUniverse;

import java.awt.Component;
import java.awt.FontMetrics;
import java.awt.Frame;
import java.awt.Window;

import javax.swing.SwingUtilities;
import javax.vecmath.Color3f;

import xmipp.ij.commons.ImagePlusLoader;
import xmipp.ij.commons.Tool;
import xmipp.ij.commons.XmippIJUtil;
import xmipp.ij.commons.XmippIJWindow;
import xmipp.ij.commons.XmippImageCanvas;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.ij.commons.XmippStackWindow;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.Param;
import xmipp.utils.XmippDialog;
import xmipp.viewer.ctf.CTFAnalyzerJFrame;
import xmipp.viewer.ctf.CTFRecalculateImageWindow;
import xmipp.viewer.ctf.TasksEngine;

/**
 * 
 * @author Juanjo Vega
 */
public class ImagesWindowFactory {

	private final static int UNIVERSE_W = 400, UNIVERSE_H = 400;

	public static void openFilesAsDefault(String filenames[], Param parameters) {
		for (int i = 0; i < filenames.length; i++) {
			openFileAsDefault(filenames[i], parameters);
		}
	}

	public static void openFileAsDefault(String filename) {
		openFileAsDefault(filename, new Param());
	}

	public static void openFileAsDefault(String filename, Param parameters) {
		try {
			if (Filename.isMetadata(filename)) {
				if (parameters.mode.equalsIgnoreCase(Param.OPENING_MODE_IMAGE))
					openFileAsImage(null, filename, parameters);
				else
					openMetadata(filename, parameters,
							Param.OPENING_MODE_GALLERY);
			} else {
				ImageGeneric img = new ImageGeneric(filename);

				if (img.isSingleImage()) {
					openFileAsImage(null, filename, parameters);
				} else if (img.isStackOrVolume()) {
					if (parameters.mode
							.equalsIgnoreCase(Param.OPENING_MODE_IMAGE))
						openFileAsImage(null, filename, parameters);
					else
						openMetadata(filename, parameters,
								Param.OPENING_MODE_GALLERY);
				}
			}
		} catch (Exception e) {
			XmippDialog.showError(null, String.format(
					"Couldn't open file: '%s'\nError: %s", filename,
					e.getMessage()));
			DEBUG.printException(e);
		}
	}

	public static void openFilesAsImages(String filenames[], Param parameters) {
		for (int i = 0; i < filenames.length; i++) {
			openFileAsImage(null, filenames[i], parameters);
		}
	}

	public static void openFileAsImage(String path) {
		openFileAsImage(null, path, new Param());
	}

	public static void openFileAsImage(Frame pframe, String filename,
			Param parameters) {
		try {
			// ImagePlus imp = openFileAsImagePlus(filename, parameters);
			ImageGeneric ig = new ImageGeneric(filename);
			ImagePlusLoader ipl = new ImagePlusLoader(ig);
			XmippIJWindow xiw = openXmippImageWindow(pframe, ipl,
					parameters.poll);
			if (parameters.mask_toolbar)
				xiw.openMaskToolbar();
		} catch (Exception e) {
			XmippDialog.showError(null, String.format(
					"Couldn't open file: '%s'\nError: %s", filename,
					e.getMessage()));
			DEBUG.printException(e);
		}
	}

	public static ImagePlus openFileAsImagePlus(String path, Param parameters)
			throws Exception {
		ImagePlus imp;
		if (Filename.isMetadata(path)) {
			MetaData md = new MetaData(path);
			imp = XmippImageConverter.readMetadataToImagePlus(
					MDLabel.MDL_IMAGE, md, parameters.useGeo, parameters.wrap);
			md.destroy();
		} else {
			imp = XmippImageConverter.loadImage(path,
					parameters.zoom > 0 ? parameters.zoom : 100);
		}
		return imp;
	}

	public static XmippIJWindow openXmippImageWindow(Window window,
			ImagePlus imp, boolean poll) {
		return openXmippImageWindow(window, new ImagePlusLoader(imp), poll);
	}

	public static XmippIJWindow openXmippImageWindow(Window window,
			ImagePlusLoader impLoader, boolean poll) {
		return openXmippImageWindow(window, impLoader, null, poll);
		
	}
	public static XmippIJWindow openXmippImageWindow(Window window,
			ImagePlusLoader impLoader, String title, boolean poll) {
		ImagePlus imp = impLoader.getImagePlus();
                
		XmippIJWindow iw;
		
		if (impLoader.isStackOrVolume())
			iw = (title != null)? new XmippStackWindow(window, impLoader, title): new XmippStackWindow(window, impLoader);
		else
			iw = (title != null )? new XmippImageWindow(impLoader, title): new XmippImageWindow(impLoader);
		SwingUtilities.invokeLater(new Worker(iw));
		return iw;
	}

	public static class Worker implements Runnable {

		XmippIJWindow iw;

		public Worker(XmippIJWindow iw) {
			this.iw = iw;
		}

		@Override
		public void run() {
			((XmippImageCanvas)  iw.getCanvas()).adjustMagnification();
			((ImageWindow) iw).setVisible(true);
		}
	}

	/**
	 * Before calling this method be sure you have constructed the proper
	 * metadata with files to be shown, mode passed will be override in
	 * parameters
	 */
	public static GalleryJFrame openMetadata(String filename, MetaData md,
			Param parameters, String mode) {

		if (parameters.mode.equalsIgnoreCase(Param.OPENING_MODE_DEFAULT))
			parameters.mode = mode;
		return new GalleryJFrame(filename, md, parameters);
	}

	public static GalleryJFrame openMetadata(String filename, Param parameters,
			String mode) throws Exception {
		return openMetadata(filename, new MetaData(filename), parameters, mode);
	}

	public static GalleryJFrame openFilesAsGallery(String filenames[],
			boolean useSameTable) throws Exception {
		return openFilesAsGallery(filenames, useSameTable, new Param());
	}

	public static GalleryJFrame openFilesAsGallery(String filenames[],
			boolean useSameTable, Param parameters) throws Exception {
		GalleryJFrame gallery = null;

		if (useSameTable) {
			MetaData md = new MetaData();
			for (int i = 0; i < filenames.length; ++i)
				md.setValueString(MDLabel.MDL_IMAGE, filenames[i],
						md.addObject());
			openMetadata(null, md, parameters, null);
		} else {
			for (int i = 0; i < filenames.length; i++) {
				gallery = openMetadata(filenames[i], parameters,
						Param.OPENING_MODE_GALLERY);
			}
		}

		return gallery;
	}

	public static void openImagePlusAs3D(ImagePlus ip) {
		try {
			// Checks if java3D is available or not.
			Class.forName("javax.media.j3d.J3DBuffer");

			new StackConverter(ip).convertToRGB();

			Image3DUniverse universe = new Image3DUniverse(UNIVERSE_W,
					UNIVERSE_H);

			// Adds the sphere image plus to universe.
			Content c = universe.addSurfacePlot(ip, new Color3f(1f, 165f / 255,
					82f / 255), "1", 50, new boolean[] { true, true, true }, 1);
			c.displayAs(Content.SURFACE);
			c.setColor(new Color3f(1f, 165f / 255, 82f / 255));

			universe.show(); // Shows...
		} catch (final ClassNotFoundException e) {
			IJ.error("Java 3D not found. Please, check your installation.");
		}
	}

	public static ImageWindow openCTFImage(ImagePlus ip, String CTFfilename,
			String PSDfilename, TasksEngine tasksEngine,
			String MicrographFilename, int row, String sortFn) {
		XmippIJUtil.showImageJ(Tool.VIEWER);// removed Toolbar.FREEROI
		return new CTFRecalculateImageWindow(ip, CTFfilename, PSDfilename,
				tasksEngine, row, sortFn);
	}

	public static void openFileAsText(String filename, Component parent) {
		TextfileJFrame frameText = new TextfileJFrame(filename);
		if (parent != null)
			frameText.setLocationRelativeTo(null);
		frameText.setVisible(true);
	}

	public static void openCTFWindow(ImagePlus imp, String CTFFilename,
			String PSDFilename) {
//		CTFProfileWindow ctfView = new CTFProfileWindow(imp, CTFFilename,
//				PSDFilename);
//		ctfView.setVisible(true);
		new CTFAnalyzerJFrame(imp, CTFFilename, PSDFilename);
	}

	public static String getSortTitle(String title, int width,
			FontMetrics fontMetrics) {
		String sort = title;
		int strlenght = fontMetrics.stringWidth(sort);
		int index = 0;

		while (strlenght > width) {
			index++;
			sort = "..." + title.substring(index);
			strlenght = fontMetrics.stringWidth(sort);
		}

		return sort;
	}
}
