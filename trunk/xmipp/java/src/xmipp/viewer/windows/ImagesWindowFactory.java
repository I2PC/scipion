/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package xmipp.viewer.windows;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageWindow;
import ij.gui.Toolbar;
import ij.io.FileInfo;
import ij.process.StackConverter;
import ij3d.Content;
import ij3d.Image3DUniverse;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Toolkit;
import java.io.File;

import javax.swing.JFrame;
import javax.swing.JRootPane;
import javax.vecmath.Color3f;

import xmipp.ij.commons.ImagePlusLoader;
import xmipp.ij.commons.Tool;
import xmipp.ij.commons.XmippIJUtil;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.ij.commons.XmippStackWindow;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.InfiniteProgressPanel;
import xmipp.utils.Param;
import xmipp.utils.XmippDialog;
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
				openMetadata(filename, parameters, Param.OPENING_MODE_GALLERY);
			} else {
				ImageGeneric img = new ImageGeneric(filename);

				if (img.isSingleImage()) {
					openFileAsImage(filename, parameters);
				} else if (img.isStackOrVolume()) {
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
			openFileAsImage(filenames[i], parameters);
		}
	}

	public static void openFileAsImage(String path) {
		openFileAsImage(path, new Param());
	}

	public static void openFileAsImage(String filename, Param parameters) {
		try {
			ImagePlus imp = openFileAsImagePlus(filename, parameters);
			openXmippImageWindow(imp, parameters.poll);
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
			imp = XmippImageConverter.readMetadataToImagePlus(MDLabel.MDL_IMAGE, md,
					parameters.useGeo, parameters.wrap);
		} else {
			imp = XmippImageConverter.loadImage(path,
					parameters.zoom > 0 ? parameters.zoom : 100);
		}
		return imp;
	}

	public static ImageWindow openXmippImageWindow(ImagePlus imp, boolean poll) {
		return openXmippImageWindow(new ImagePlusLoader(imp), poll);
	}
	
	public static ImageWindow openXmippImageWindow(ImagePlusLoader impLoader, boolean poll) {
		ImagePlus imp = impLoader.getImagePlus();
		ImageWindow iw;
		if (imp.getStackSize() > 1)
			iw = new XmippStackWindow(impLoader);
		else
			iw = new XmippImageWindow(impLoader);
		iw.setVisible(true);
		return iw;
	}

	/**
	 * Before calling this method be sure you have constructed the proper
	 * metadata with files to be shown, mode passed will be override in
	 * parameters
	 */
	public static JFrameGallery openMetadata(String filename, MetaData md,
			Param parameters, String mode) {

		if (parameters.mode.equalsIgnoreCase(Param.OPENING_MODE_DEFAULT))
			parameters.mode = mode;
		return new JFrameGallery(filename, md, parameters);
	}

	public static JFrameGallery openMetadata(String filename, Param parameters,
			String mode) throws Exception {
		return openMetadata(filename, new MetaData(filename), parameters, mode);
	}

	public static JFrameGallery openFilesAsGallery(String filenames[],
			boolean useSameTable) throws Exception {
		return openFilesAsGallery(filenames, useSameTable, new Param());
	}

	public static JFrameGallery openFilesAsGallery(String filenames[],
			boolean useSameTable, Param parameters) throws Exception {
		JFrameGallery gallery = null;

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
			String MicrographFilename, int row) {
		XmippIJUtil.showImageJ(Tool.VIEWER);// removed Toolbar.FREEROI
		return new CTFRecalculateImageWindow(ip, CTFfilename, PSDfilename,
				tasksEngine, row);
	}

	public static void openFileAsText(String filename, Component parent) {
		JFrameTextfile frameText = new JFrameTextfile(filename);
		if (parent != null)
			frameText.setLocationRelativeTo(null);
		frameText.setVisible(true);
	}

	public static void openCTFWindow(ImagePlus imp, String CTFFilename,
			String PSDFilename) {
		CTFProfileWindow ctfView = new CTFProfileWindow(imp, CTFFilename,
				PSDFilename);
		ctfView.setVisible(true);
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
