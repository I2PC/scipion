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
 *   - Why?
 * Centralize data file operations (read, write...) in a single module (class)
 * One couple of methods per file type
 * 
 *   - Alternatives
 * Make specialized subclasses, one per file type. Possible advantages: 
 * > smaller source files 
 */
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;

public class TiltSeriesIO {
	// this class uses JNI mechanisms to call C++ code
	// TODO: file formats: spider (MRC and XMP already work)
	static {
		// loads the XMIPP C++ library for data
		System.loadLibrary("XmippDataJava");
	}
	
	// if an image is bigger than this threshold, resize it to this size
	// This optimization makes sense since for interactive use no more detail is needed
	public static Dimension resizeThreshold = new Dimension(400,400);

	/**
	 * @param model
	 *            stores relevant data following MVC paradigm
	 * @throws java.io.IOException
	 * @throws InterruptedException
	 *             Threads issues
	 */
	public static void read(TomoData model) throws java.io.IOException,
			InterruptedException {

		String path = model.getFilePath();

		if ((path == null) || (path.equals("")))
			throw new IOException("Empty path");

		// path=path.toLowerCase();

		// right now, identify file type by its extension
		// should check also uppercase like .SEL
		if (path.endsWith(".mrc") || path.endsWith(".mrcs")) {
			readImages(path, model);
		} else if (path.endsWith(".sel")) {
			readSel(path, model);
		} else
			throw new IOException("Unknown file type");

	}

	public static void write(TomoData model) {
		// TODO: write(model) - use C++ version
		// TODO: write(model) - check destination path, so we write to the desired file (and not to some unexpected other)
		// TODO: write(model) - fill the ENABLED field in sel file (so discarded projections are marked as disabled)
		String path = model.getFilePath();

		if (path.endsWith(".mrc")) {
			MrcIOJava.writeMRC(model);
		}
	}

	/**
	 * @param path
	 *            absolute path to selfile
	 * @param model
	 * @throws IOException
	 */
	private static void readSel(String path, TomoData model) throws IOException {
		// TODO: readSel - read tilt angles (from Metadata)
		// TODO: readSel - CRASH - old selfiles use relative paths, new selfiles use absolute paths
		
		// files inside selfiles are relative to selfile location
		// baseDir only needed with ImageJ FileDialog
		// XmippTomo FileDialog returns absolute path
		//String baseDir = new File(path).getParent();
		String baseDir=null;

		MetaData md = new MetaData();

		FileName fn = new FileName(path);
		// Xmipp_Tomo.debug(fn.getBaseName());
		md.read(fn);
		md.iteratorBegin();
		while (!md.iteratorEnd()) {
			// store in filenameString the field MDL_IMAGE of each row of the
			// selfile
			String filenameString[] = new String[1];
			md.getStrFromValue(MDLabel.MDL_IMAGE, filenameString);
			try {
				// Xmipp_Tomo.debug(filenameString[0]);
				ImagePlusC ip = readImage(filenameString[0], baseDir, false);
				postReadImage(model, ip);
			} catch (IOException ex) {
				// don't throw exception, so it tries to read next image
				Xmipp_Tomo.debug("readimages - ", ex);
			}
			md.iteratorNext();
		}

		postReadImages(model);
	}

	private static ImagePlusC readImage(String path, String base, boolean headerOnly)
			throws IOException {
		return readImage(path, base, headerOnly, -1);
	}

	/**
	 * Read an image using the C++ Xmipp code (see java_wrapper.h)
	 * 
	 * @param path
	 *            to file - if relative, base should not be null
	 * @param base
	 *            base directory (for relative paths)
	 * @param headerOnly
	 *            - read only the header of the image (vs whole image data)
	 * @param slice
	 *            - slice number (0..N-1)
	 * @return image as ImagePlusC
	 * @throws IOException
	 */
	private static ImagePlusC readImage(String path, String base, boolean headerOnly,
			int slice) throws IOException {
		String filepath = path;
		if (base != null)
			filepath = base + "/" + path;

		// Xmipp_Tomo.debug(filepath);

		ImagePlusC ipc = new ImagePlusC();
		ipc.setFilename(filepath);
		ipc.setReadHeaderOnly(headerOnly);
		ipc.setSlice(slice);
		int err = XmippData.readImage(ipc);

		if (err != 0)
			throw new IOException("Image.read error " + err);

		return ipc;
	}

	/**
	 * Actions required after reading 1 ImagePlusC
	 * 
	 * @param model where the postprocessed data is saved
	 * @param img
	 * @depends readImage
	 */
	private static void postReadImage(TomoData model, ImagePlusC img) {
		ImagePlus image = convert(img);

		model.setOriginalHeight(img.getHeight());
		model.setOriginalWidth(img.getWidth());
		model.updateMinMax(image);

		// resize/scale - use an aux image processor for all projections
		ImageProcessor ipresized = null, ip = image.getProcessor();

		if (shouldResize(img.getWidth(),img.getHeight())) {
			ipresized = ip.resize(resizeThreshold.width,resizeThreshold.height);
			model.addProjection(ipresized);
			model.setResized(true);
		} else
			model.addProjection(ip);
	}
	
	/**
	 * Actions required after reading all the images
	 * 
	 * @param model
	 */
	private static void postReadImages(TomoData model) {
		model.lastImageLoaded();
	}

	
	/**
	 * @param img
	 * @return
	 */
	private static ImagePlus convert(ImagePlusC img) {

		// Creates image
		FloatProcessor convertedImageProcessor = new FloatProcessor(img.getWidth(), img.getHeight());

		if (img.getReadHeaderOnly() == false)
			for (int x = 0; x < img.getWidth(); x++)
				for (int y = 0; y < img.getHeight(); y++)
					convertedImageProcessor.setf(x, y, (float) img.getPixel(x, y));

		// normalize the image - done by default when calling FloatProcessor(array), here is our responsibility
		convertedImageProcessor.findMinAndMax();
		
		ImagePlus imagePlus = new ImagePlus("ImageDouble", convertedImageProcessor);
		// imagePlus.getProcessor().setPixels(img.getImage());
		imagePlus.setDimensions(0, img.getNImages(), 0);
		return imagePlus;
	}

	/**
	 * @param path
	 *            absolute path to file
	 * @param model
	 * @throws java.io.IOException
	 * @throws InterruptedException
	 * @throws OutOfMemoryError
	 */
	private static void readImages(String path, TomoData model)
			throws java.io.IOException, InterruptedException, OutOfMemoryError {
		// files are relative to selfile location
		ImagePlusC img = readImage(path, null, true);
		int nImages = img.getNImages();

		model.setHeight(img.getHeight());
		model.setWidth(img.getWidth());

		// load images into model
		// Xmipp_Tomo.debug("Number of images: " + nImages);
		for (int i = 1; i <= nImages; i++) {
			if (model.isLoadCanceled())
				break;
			try {
				ImagePlusC image = readImage(path, null, false, i - 1);
				postReadImage(model, image);
			} catch (IOException ex) {
				// don't throw exception, so it tries to read next image
				Xmipp_Tomo.debug("readimages - ", ex);
			}
			// uncomment for concurrency tests
			// if(Xmipp_Tomo.TESTING != 0)
			// Thread.sleep(250);
		}
		postReadImages(model);

		// read tilt angles
		String tltFilePath = getTiltFilePath(path);
		try {
			readTiltAngles(tltFilePath, model);
		} catch (FileNotFoundException ex) {
			// model.emptyTiltAngles(model.getNumberOfProjections());
			Xmipp_Tomo.debug("readImages - problem reading tlt. "
					+ ex.toString());
		}

	}



	/**
	 * @param tltFilePath
	 * @param model
	 * @throws java.io.IOException
	 */
	public static void readTiltAngles(String tltFilePath, TomoData model)
			throws java.io.IOException {
		// .tlt syntax: one float angle per line, stored as text
		
		BufferedReader brin = new BufferedReader(new FileReader(tltFilePath));
		String line = null;
		while ((line = brin.readLine()) != null) {
			model.addTiltAngle(new Float(line));
		}
	}

	/**
	 * TODO: writeTiltAngles - test
	 * @param tltFilePath
	 * @param model
	 * @throws java.io.IOException
	 */
	public static void writeTiltAngles(String tltFilePath, TomoData model)
			throws IllegalArgumentException {
		// .tlt syntax: one float angle per line, stored as text
		if ((tltFilePath == null) || ("".equals(tltFilePath))) {
			throw new IllegalArgumentException("Null path");
		}

		File f = new File(tltFilePath);

		/*
		 * Extra checks for files that already exist, and hence are going to be
		 * "rewritten" if (!f.isFile()) { throw new
		 * IllegalArgumentException("Path is a directory: " + tltFilePath); } if
		 * (!f.canWrite()) { throw new
		 * IllegalArgumentException("File cannot be written: " + tltFilePath); }
		 */

		// Writer output = new BufferedWriter(new FileWriter(f));
		PrintWriter output = null;
		try {
			output = new PrintWriter(f);
			Iterator<Float> i = model.getTiltAnglesIterator();
			do {
				output.println(i.next());
			} while (i.hasNext());

		} catch (FileNotFoundException ex) {
			Xmipp_Tomo.debug("writeTlt - file not found");
		} finally {
			output.close();
		}

	}



	protected static String getTiltFilePath(String path) {
		return changeExtension(path, ".tlt");
	}

	// Copyleft 2003 Fred Swartz
	private static String changeExtension(String originalName,
			String newExtension) {
		int lastDot = originalName.lastIndexOf(".");
		if (lastDot != -1) {
			return originalName.substring(0, lastDot) + newExtension;
		} else {
			return originalName + newExtension;
		}
	}


	public static boolean shouldResize(int width, int height){
		// TODO: shouldResize - when applying changes to the original image, don't resize
		return (width > resizeThreshold.width ) || (height > resizeThreshold.height);
	}


}
