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
import xmipp.ij.commons.XmippIJWindow;
import xmipp.ij.commons.XmippImageCanvas;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.ij.commons.XmippMenuBar;
import xmipp.ij.commons.XmippStackWindow;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.Params;
import xmipp.utils.ScipionParams;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippMessage;
import xmipp.viewer.scipion.ScipionGalleryData;
import xmipp.viewer.scipion.ScipionGalleryJFrame;
import xmipp.viewer.scipion.ScipionMetaData;

/**
 * 
 * @author Juanjo Vega
 */
public class ImagesWindowFactory {

	private final static int UNIVERSE_W = 400, UNIVERSE_H = 400;

	public static void openFilesAsDefault(String filenames[], Params parameters) {
		for (int i = 0; i < filenames.length; i++) {
			openFileAsDefault(filenames[i], parameters);
		}
	}

	public static void openFileAsDefault(String filename) {
		
		Params params = new Params();
		
		if (ScipionMetaData.isScipionMetaData(filename)) {
			params = new ScipionParams();
		}
		
		openFileAsDefault(filename, params);
	}

	public static void openFileAsDefault(String filename, Params parameters) {
		try {
			if (Filename.isMetadata(filename)) {
				if (parameters.mode.equalsIgnoreCase(Params.OPENING_MODE_IMAGE))
					openFileAsImage(null, filename, parameters);
				else
					openMetadata(filename, parameters,
							Params.OPENING_MODE_GALLERY);
			} else {
				ImageGeneric img = new ImageGeneric(filename);

				if (img.isSingleImage()) {
					openFileAsImage(null, filename, parameters);
				} else if (img.isStackOrVolume()) {
					if (parameters.mode
							.equalsIgnoreCase(Params.OPENING_MODE_IMAGE))
						openFileAsImage(null, filename, parameters);
					else
						openMetadata(filename, parameters, Params.OPENING_MODE_GALLERY);
				}
			}
		} catch (Exception e) {
			XmippDialog.showError(null, String.format(
					"Couldn't open file: '%s'\nError: %s", filename,
					e.getMessage()));
			e.printStackTrace();
		}
	}

	public static void openFilesAsImages(String filenames[], Params parameters) {
		for (int i = 0; i < filenames.length; i++) {
			openFileAsImage(null, filenames[i], parameters);
		}
	}

	public static void openFileAsImage(String path) {
		openFileAsImage(null, path, new Params());
	}

	public static void openFileAsImage(Frame pframe, String filename,
			Params parameters) {
		try {
			ImagePlusLoader ipl = new ImagePlusLoader(filename, false);
			XmippIJWindow xiw = openXmippImageWindow(ipl,
					parameters);
			if (parameters.mask_toolbar)
				xiw.openMaskToolbar();
		} catch (Exception e) {
			XmippDialog.showError(null, String.format(
					"Couldn't open file: '%s'\nError: %s", filename,
					e.getMessage()));
			e.printStackTrace();
		}
	}

	public static ImagePlus openFileAsImagePlus(String path, Params parameters)
			throws Exception {
		ImagePlus imp;
		if (Filename.isMetadata(path)) {
			MetaData md = new MetaData(path);
			imp = XmippImageConverter.readMetadataToImagePlus(
					MDLabel.MDL_IMAGE, md, parameters.useGeo, parameters.wrap, parameters.inverty);
			md.destroy();
		} else {
			imp = XmippImageConverter.loadImage(path,
					parameters.zoom != null ? parameters.zoom : 100);
		}
		return imp;
	}

	public static XmippIJWindow openXmippImageWindow(
			ImagePlus imp, Params params) {
		return openXmippImageWindow(new ImagePlusLoader(imp), params);
	}

	public static XmippIJWindow openXmippImageWindow(
			ImagePlusLoader impLoader, Params parameters) {
		return openXmippImageWindow(impLoader, impLoader.getName(), parameters);
		
	}
	public static XmippIJWindow openXmippImageWindow(
			ImagePlusLoader impLoader, String title, Params parameters) {
                
		XmippIJWindow iw;
		
		if (impLoader.isStackOrVolume())
			iw = (title != null)? new XmippStackWindow(impLoader, title, parameters): new XmippStackWindow(impLoader, parameters);
		else
			iw = (title != null )? new XmippImageWindow(impLoader, title, parameters): new XmippImageWindow(impLoader, parameters);
                
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
                        Frame frame = (ImageWindow) iw;
			frame.setVisible(true);
                        
		}
	}

	/**
	 * Before calling this method be sure you have constructed the proper
	 * metadata with files to be shown, mode passed will be override in
	 * parameters
	 */
	public static GalleryJFrame openMetadata(String filename, MetaData md,
			Params parameters, String mode) {

		if (parameters.mode.equalsIgnoreCase(Params.OPENING_MODE_DEFAULT))
			parameters.mode = mode;
		return new GalleryJFrame(md, parameters);
	}

	public static void openMetadata(String filename, Params parameters,
			String mode) throws Exception 
        {
                if(filename.endsWith(".sqlite") || filename.endsWith(".db"))
                    openScipionMetadata(filename, parameters);
                else
                    openMetadata(filename, new MetaData(filename), parameters, mode);
	}
        
        public static void openScipionMetadata(String filename, Params parameters) throws Exception 
        {
                if(!(filename.endsWith(".sqlite") || filename.endsWith(".db")))
                    throw new IllegalArgumentException(XmippMessage.getIllegalValueMsg("scipion metadata", filename));
                else
                {
                    ScipionMetaData md;
                    if(filename.contains("@"))
                    {
                        int sep = filename.indexOf('@');
                        final String preffix = filename.substring(0, sep);
                        filename = filename.substring(sep + 1);
                        md = new ScipionMetaData(filename);
                        final ScipionGalleryData data = new ScipionGalleryData(null, parameters, md);
                        SwingUtilities.invokeLater(new Runnable() {

                            @Override
                            public void run() {
                                ScipionGalleryJFrame frame = new ScipionGalleryJFrame(data);
                                frame.selectBlock(preffix);
                            }
                        });
                        
                    }
                    else
                    {
                        md = new ScipionMetaData(filename);
                        ScipionGalleryData data = new ScipionGalleryData(null, parameters, md);
                        new ScipionGalleryJFrame(data);
                    }
                }
	}
        
       
        

	public static void openFilesAsGallery(String filenames[],
			boolean useSameTable) throws Exception {
		openFilesAsGallery(filenames, useSameTable, new Params());
	}

	public static void openFilesAsGallery(String filenames[],
			boolean useSameTable, Params parameters) throws Exception {
		GalleryJFrame gallery = null;

		if (useSameTable) {
			MetaData md = new MetaData();
			for (int i = 0; i < filenames.length; ++i)
				md.setValueString(MDLabel.MDL_IMAGE, filenames[i],
						md.addObject());
			openMetadata(null, md, parameters, null);
		} else {
			for (int i = 0; i < filenames.length; i++) {
				openMetadata(filenames[i], parameters,
						Params.OPENING_MODE_GALLERY);
			}
		}

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

	

	public static void openFileAsText(String filename, Component parent) {
		TextfileJFrame frameText = new TextfileJFrame(filename);
		if (parent != null)
			frameText.setLocationRelativeTo(null);
		frameText.setVisible(true);
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
