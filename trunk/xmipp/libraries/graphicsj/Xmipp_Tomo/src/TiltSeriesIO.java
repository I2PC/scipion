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
import ij.ImageStack;
import ij.io.FileInfo;
import ij.io.FileOpener;
import ij.io.ImageReader;
import ij.io.ImageWriter;
import ij.measure.Calibration;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.image.ColorModel;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.io.Writer;
import java.util.Iterator;
import java.util.Properties;

/**
 * Collection of methods for handling I/O - make all methods static?
 */
public class TiltSeriesIO {
	static {
		// loads "XmippDataJava"
		System.loadLibrary("XmippDataJava");
	}

	public static int MRC_HEADER_SIZE = 1024;
	private boolean resize = true;

	// Header of a data file
	private class ByteHeader {
		private byte[] bytes;
		private boolean littleEndian = false;

		public ByteHeader(int size) {
			bytes = new byte[size];
		}

		// fill bytes buffer with data from path
		public void read(String path) throws java.io.IOException {
			RandomAccessFile f = new RandomAccessFile(path, "r");

			for (int i = 0; i < bytes.length; i++) {
				bytes[i] = f.readByte();
			}
			f.close();
		}

		public void write(DataOutputStream out) throws java.io.IOException {
			for (byte b : bytes) {
				out.writeByte(b);
			}
		}

		/**
		 * @param index
		 *            - position in the header, structured as if it were a
		 *            series of integers (that is, the real position in the byte
		 *            header will be the index times the size of an integer)
		 * @param value
		 */
		public void setInt(int index, int value) {
			int i = index * 4;
			byte b1 = (byte) (value >> 24);
			byte b2 = (byte) ((value << 8) >> 24);
			byte b3 = (byte) ((value << 16) >> 24);
			byte b4 = (byte) ((value << 24) >> 24);
			if (isLittleEndian()) {
				bytes[i] = b4;
				bytes[i + 1] = b3;
				bytes[i + 2] = b2;
				bytes[i + 3] = b1;
			} else {
				bytes[i] = b1;
				bytes[i + 1] = b2;
				bytes[i + 2] = b3;
				bytes[i + 3] = b4;
			}
		}

		// return 4 bytes at position index of the buffer (interpreted as an
		// integer series)
		public int getInt(int index) {
			byte b1 = bytes[index];
			byte b2 = bytes[index + 1];
			byte b3 = bytes[index + 2];
			byte b4 = bytes[index + 3];
			if (isLittleEndian()) {
				return ((((b4 & 0xff) << 24) | ((b3 & 0xff) << 16)
						| ((b2 & 0xff) << 8) | (b1 & 0xff)));
			}
			return ((((b1 & 0xff) << 24) | ((b2 & 0xff) << 16)
					| ((b3 & 0xff) << 8) | (b4 & 0xff)));
		}

		/**
		 * @return is littleEndian?
		 */
		public boolean isLittleEndian() {
			return littleEndian;
		}

		/**
		 * @param littleEndian
		 *            the littleEndian to set
		 */
		public void setLittleEndian(boolean littleEndian) {
			this.littleEndian = littleEndian;
		}
	}

	public TiltSeriesIO(boolean resize) {
		this.resize = resize;
	}

	/**
	 * @param path
	 *            file to read
	 * @param model
	 *            stores relevant data following MVC paradigm
	 * @throws java.io.IOException
	 * @throws InterruptedException
	 *             Threads issues
	 */
	public void read(TomoData model) throws java.io.IOException,
			InterruptedException {

		String path = model.getFilePath();

		if ((path == null) || (path.equals("")))
			throw new IOException("Empty path");

		// path=path.toLowerCase();

		// right now, identify file type by its extension
		if (path.endsWith(".mrc") || path.endsWith(".mrcs")) {
			readImages(path, model);

		} else if (path.endsWith(".sel")) {
			readSel(path, model);
		}else
			throw new IOException("Unknown file type");
		/*
		 * else if (path.endsWith(".mrcs")) { readMRCS(path,model); }
		 */
		/*
		 * else if (path.toLowerCase().endsWith(".spi") ||
		 * path.toLowerCase().endsWith(".xmp")||
		 * path.toLowerCase().endsWith(".vol")) { imp = (ImagePlus)
		 * IJ.runPlugIn("Spider_Reader", path); if (imp == null) { //width =
		 * PLUGIN_NOT_FOUND; } if (imp != null && imp.getWidth() == 0) { imp =
		 * null; } }else if (path.endsWith(".sel")) { imp = (ImagePlus)
		 * IJ.runPlugIn("Sel_Reader", path); if (imp == null) { //width =
		 * PLUGIN_NOT_FOUND; } if (imp != null && imp.getWidth() == 0) { imp =
		 * null; } }
		 */

	}

	public void write(TomoData model) {
		String path = model.getFilePath();

		if (path.endsWith(".mrc")) {
			writeMRC(model);
		}
	}

	/**
	 * selfile includes tilt angles
	 * 
	 * @param path
	 *            absolute path to selfile
	 * @param model
	 * @throws IOException
	 */
	private void readSel(String path, TomoData model) throws IOException {

		// files inside selfiles are relative to selfile location
		String baseDir = new File(path).getParent();

		MetaData md = new MetaData();

		FileName fn = new FileName(path);
		Xmipp_Tomo.debug(fn.getBaseName());
		md.read(fn);
		int k = md.iteratorBegin();
		while (!md.iteratorEnd()) {
			// store in filenameString the field MDL_IMAGE of each row of the
			// selfile
			String filenameString[] = new String[1];
			md.getStrFromValue(MDLabel.MDL_IMAGE, filenameString);

			ImagePlusC ip = readImage(filenameString[0], baseDir, false);
			postReadImage(model, ip);
			k = md.iteratorNext();
		}

		postReadImages(model);
	}

	private ImagePlusC readImage(String path, String base, boolean headerOnly)
			throws IOException {
		return readImage(path, base, headerOnly, -1);
	}

	private ImagePlusC readImage(String path, String base, boolean headerOnly,
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

	private void postReadImage(TomoData model, ImagePlusC img) {
		ImagePlus image=convert(img);
		// resize/scale - use an aux image processor for all projections
		ImageProcessor ipresized = null, ip = image.getProcessor();

		if (isResize()) {
			ipresized = ip.resize(Xmipp_Tomo.resizeThreshold.width,
					Xmipp_Tomo.resizeThreshold.height);
			model.addProjection(ipresized);
		} else
			model.addProjection(ip);
	}

	private void postReadImages(TomoData model) {
		if (isResize()) {
			// Adjust file info to scaled size
			//fi.width = Xmipp_Tomo.resizeThreshold.width;
			//fi.height = Xmipp_Tomo.resizeThreshold.height;
			model.setResized(true);
		}

		model.lastImageLoaded();
	}


	private ImagePlus convert(ImagePlusC img) {

		// Creates image
		FloatProcessor ip = new FloatProcessor(img.getWidth(), img.getHeight());
		
		if(img.getReadHeaderOnly() == false)
			for (int x = 0; x < img.getWidth(); x++)
				for (int y = 0; y < img.getHeight(); y++)
					ip.setf(x, y, (float) img.getPixel(x, y));
		
		ImagePlus imagePlus = new ImagePlus("ImageDouble", ip);
		// imagePlus.getProcessor().setPixels(img.getImage());
		imagePlus.setDimensions(0, img.getNImages(), 0);
		return imagePlus;

	}

	/**
	 * Read images from a file like MRC or SPI
	 * 
	 * @param path
	 *            absolute path to file
	 * @param model
	 * @throws java.io.IOException
	 * @throws InterruptedException
	 * @throws OutOfMemoryError
	 */
	private void readImages(String path, TomoData model)
			throws java.io.IOException, InterruptedException, OutOfMemoryError {
		// files are relative to selfile location
		ImagePlusC img = readImage(path, null, true);
		int nImages = img.getNImages();

		// Save original size to model
		model.setHeight(img.getHeight());
		model.setWidth(img.getWidth());

		// load images
		Xmipp_Tomo.debug("Number of images: " + nImages);
		for (int i = 1; i <= nImages; i++) {
			if (model.isLoadCanceled())
				break;

			ImagePlusC ip = readImage(path, null, false, i - 1);
			postReadImage(model, ip);

			// uncomment for concurrency tests
			// if(Xmipp_Tomo.TESTING != 0)
			// Thread.sleep(250);
		}
		postReadImages(model);
		
		// read tilt angles
		String tltFilePath = getTltPath(path);
		try {
			readTlt(tltFilePath, model);
		} catch (FileNotFoundException ex) {
			// model.emptyTiltAngles(model.getNumberOfProjections());
		}
		
	}

	private void readMRCJava(String path, TomoData model)
			throws java.io.IOException, InterruptedException, OutOfMemoryError {

		// Parse header
		ByteHeader header = new ByteHeader(MRC_HEADER_SIZE);
		header.read(path);
		FileInfo fi = prepareFileInfo(path, header);

		// read each slice, resize/scale and add it to the stack
		FileOpener fo = new FileOpener(fi);

		ColorModel cm = fo.createColorModel(fi);

		// ImageStack stack = new ImageStack(fi.width, fi.height, cm);
		long skip = fi.longOffset > 0 ? fi.longOffset : fi.offset;
		Object pixels;

		ImageReader reader = new ImageReader(fi);
		InputStream is = fo.createInputStream(fi);
		if (is == null)
			throw new IOException(
					"TiltSeriesOpener.readMRC - Null input stream");

		// Save original size to model
		model.setHeight(fi.height);
		model.setWidth(fi.width);

		// load images
		Xmipp_Tomo.debug("Number of images: " + fi.nImages);
		for (int i = 1; i <= fi.nImages; i++) {
			if (model.isLoadCanceled())
				break;
			pixels = reader.readPixels(is, skip);
			if (pixels == null)
				break;

			// resize/scale - use an aux image processor for all projections
			ImageProcessor ip = null, ipresized = null;

			// create the proper(byte/short/float...) ImageProcessor initialized
			// with the pixels
			switch (fi.fileType) {
			case FileInfo.GRAY8:
				ip = new ByteProcessor(fi.width, fi.height, (byte[]) pixels, cm);
				break;
			case FileInfo.GRAY16_UNSIGNED:
				ip = new ShortProcessor(fi.width, fi.height, (short[]) pixels,
						cm);
				break;
			case FileInfo.GRAY32_FLOAT:
				ip = new FloatProcessor(fi.width, fi.height, (float[]) pixels,
						cm);
				break;
			}

			if (isResize())
				ipresized = ip.resize(Xmipp_Tomo.resizeThreshold.width,
						Xmipp_Tomo.resizeThreshold.height);

			// stack.addSlice(null, ip);

			if (isResize())
				model.addProjection(ipresized);
			else
				model.addProjection(ip);

			skip = fi.gapBetweenImages;

			// uncomment for concurrency tests
			// if(Xmipp_Tomo.TESTING != 0)
			// Thread.sleep(250);
		}
		is.close();

		if (fi.info != null)
			model.getImage().setProperty("Info", fi.info);

		if (isResize()) {
			// Adjust file info to scaled size
			fi.width = Xmipp_Tomo.resizeThreshold.width;
			fi.height = Xmipp_Tomo.resizeThreshold.height;
			model.setResized(true);
		}

		model.getImage().setFileInfo(fi);

		/*
		 * if (stack.getSize()==0) throw new IOException("Empty stack");
		 */

		/*
		 * if (fi.sliceLabels!=null && fi.sliceLabels.length<=stack.getSize()) {
		 * for (int i=0; i<fi.sliceLabels.length; i++)
		 * stack.setSliceLabel(fi.sliceLabels[i], i+1); }
		 */

		setCalibration(model.getImage(), fi, fo);
		ImageProcessor ip = model.getImage().getProcessor();
		if (ip.getMin() == ip.getMax()) // find stack min and max if first slice
										// is blank
			setStackDisplayRange(model.getImage());

		model.lastImageLoaded();

		// read tilt angles
		String tltFilePath = getTltPath(path);
		try {
			readTlt(tltFilePath, model);
		} catch (FileNotFoundException ex) {
			// model.emptyTiltAngles(model.getNumberOfProjections());
		}
	}

	private void writeMRC(TomoData model) {
		boolean littleEndian = true;
		ByteHeader header = createMRCHeader(model, littleEndian);
		ImagePlus img = model.getImage();
		FileInfo fi = img.getFileInfo();
		fi.fileFormat = fi.RAW;
		fi.intelByteOrder = littleEndian;
		fi.fileName = model.getFileName();
		fi.directory = model.getDirectory();
		ImageWriter file = new ImageWriter(fi);
		try {
			DataOutputStream out = new DataOutputStream(
					new BufferedOutputStream(new FileOutputStream(model
							.getFilePath())));
			header.write(out);
			file.write(out);
			out.close();
		} catch (IOException ioe) {
			Xmipp_Tomo.debug("TiltSeriesOpener.writeMRC" + ioe);
		}

		// write tilt file
		String tltFilePath = getTltPath(model.getFilePath());
		writeTlt(tltFilePath, model);
	}

	private ByteHeader createMRCHeader(TomoData model, boolean littleEndian) {
		ImagePlus img = model.getImage();
		ImageProcessor ip = img.getProcessor();

		int mode = 0;
		if (ip instanceof ByteProcessor) {
			mode = 0;
		} else if (ip instanceof ShortProcessor) {
			mode = 1;
		} else if (ip instanceof FloatProcessor) {
			mode = 2;
		}

		double sum = 0;
		int nbPix = 0;
		float min = Float.NaN;
		float max = Float.NaN;
		float value;
		int mini = 0;
		int maxi = 0;
		int avgi = 0;
		for (int z = 0; z < img.getNSlices(); z++) {
			ip = img.getImageStack().getProcessor(z + 1);
			for (int y = 0; y < img.getHeight(); y++) {
				for (int x = 0; x < img.getWidth(); x++) {
					value = ip.getPixelValue(x, y);
					sum += value;
					nbPix++;
					if (Float.isNaN(min)) {
						min = value;
					}
					if (Float.isNaN(max)) {
						max = value;
					}
					if (max < value) {
						max = value;
					} else if (min > value) {
						min = value;
					}
				}
			}
		}
		float avg = (float) (sum / nbPix);
		maxi = Float.floatToRawIntBits(max);
		mini = Float.floatToRawIntBits(min);
		avgi = Float.floatToRawIntBits(avg);

		ByteHeader header = new ByteHeader(1024);
		header.setLittleEndian(littleEndian);
		header.setInt(0, img.getWidth());
		header.setInt(1, img.getHeight());
		header.setInt(2, img.getNSlices());
		header.setInt(3, mode);
		header.setInt(16, 1);
		header.setInt(17, 2);
		header.setInt(18, 3);
		header.setInt(19, mini);
		header.setInt(20, maxi);
		header.setInt(21, avgi);
		return header;

	}

	/**
	 * get MRC DM3 specific information from header into FileInfo
	 * 
	 * @param path
	 * @param header
	 *            of file "path"
	 * @return FileInfo with data from header
	 * @throws java.io.IOException
	 */
	public FileInfo prepareFileInfo(String path, ByteHeader header)
			throws java.io.IOException {
		// set filename
		FileInfo fi = new FileInfo();
		fi.fileFormat = FileInfo.RAW;
		File dest = new File(path);
		fi.fileName = dest.getName();
		fi.directory = dest.getParent();

		int mode = header.getInt(3 * 4);
		if (mode == 0) {
			fi.fileType = FileInfo.GRAY8;
			fi.width = header.getInt(0);
			fi.height = header.getInt(1 * 4);
			fi.nImages = header.getInt(2 * 4);
			if (fi.width < 0 || fi.height < 0 || fi.nImages < 0) {
				header.setLittleEndian(true);
			}
		} else if (mode == 1) {
			fi.fileType = FileInfo.GRAY16_UNSIGNED;
		} else if (mode == 2) {
			fi.fileType = FileInfo.GRAY32_FLOAT;
		} else {
			// it is not a big endian data
			header.setLittleEndian(true);
			mode = header.getInt(3 * 4);
			if (mode == 1) {
				fi.fileType = FileInfo.GRAY16_UNSIGNED;
			} else if (mode == 2) {
				fi.fileType = FileInfo.GRAY32_FLOAT;
			} else {
				throw new IOException("Unimplemented ImageData mode=" + mode
						+ " in MRC file.");
			}
		}

		fi.intelByteOrder = header.isLittleEndian();
		// duplicated when mode = 0 ?
		fi.width = header.getInt(0);
		fi.height = header.getInt(1 * 4);
		fi.nImages = header.getInt(2 * 4);

		int extrabytes = header.getInt(23 * 4);
		fi.offset = 1024 + extrabytes;
		/*
		 * if (extrabytes == 131072) { IJ.write("FEI");
		 * 
		 * FEI = true; }
		 */
		return fi;
	}

	/**
	 * .tlt syntax: one float angle per line, stored as text
	 * 
	 * @param tltFilePath
	 * @param model
	 * @throws java.io.IOException
	 */
	public void readTlt(String tltFilePath, TomoData model)
			throws java.io.IOException {

		BufferedReader brin = new BufferedReader(new FileReader(tltFilePath));
		String line = null;
		while ((line = brin.readLine()) != null) {
			model.addTiltAngle(new Float(line));
		}
	}

	/**
	 * .tlt syntax: one float angle per line, stored as text
	 * 
	 * @param tltFilePath
	 * @param model
	 * @throws java.io.IOException
	 */
	public void writeTlt(String tltFilePath, TomoData model)
			throws IllegalArgumentException {
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

	// code from FileOpener - due to access modifiers restrictions...
	private void setCalibration(ImagePlus imp, FileInfo fi, FileOpener fo) {

		if (fi.fileType == FileInfo.GRAY16_SIGNED) {
			if (IJ.debugMode)
				IJ.log("16-bit signed");
			double[] coeff = new double[2];
			coeff[0] = -32768.0;
			coeff[1] = 1.0;
			imp.getLocalCalibration().setFunction(Calibration.STRAIGHT_LINE,
					coeff, "gray value");
		}

		Properties props = fo.decodeDescriptionString(fi);
		Calibration cal = imp.getCalibration();
		boolean calibrated = false;
		if (fi.pixelWidth > 0.0 && fi.unit != null) {
			cal.pixelWidth = fi.pixelWidth;
			cal.pixelHeight = fi.pixelHeight;
			cal.pixelDepth = fi.pixelDepth;
			cal.setUnit(fi.unit);
			calibrated = true;
		}

		if (fi.valueUnit != null) {
			int f = fi.calibrationFunction;
			if ((f >= Calibration.STRAIGHT_LINE && f <= Calibration.RODBARD2 && fi.coefficients != null)
					|| f == Calibration.UNCALIBRATED_OD) {
				boolean zeroClip = props != null
						&& props.getProperty("zeroclip", "false")
								.equals("true");
				cal.setFunction(f, fi.coefficients, fi.valueUnit, zeroClip);
				calibrated = true;
			}
		}

		/*
		 * if (calibrated) checkForCalibrationConflict(imp, cal);
		 */

		if (fi.frameInterval != 0.0)
			cal.frameInterval = fi.frameInterval;

		if (props == null)
			return;

		cal.xOrigin = getDouble(props, "xorigin");
		cal.yOrigin = getDouble(props, "yorigin");
		cal.zOrigin = getDouble(props, "zorigin");
		cal.info = props.getProperty("info");

		cal.fps = getDouble(props, "fps");
		cal.loop = getBoolean(props, "loop");
		cal.frameInterval = getDouble(props, "finterval");
		cal.setTimeUnit(props.getProperty("tunit", "sec"));

		double displayMin = getDouble(props, "min");
		double displayMax = getDouble(props, "max");
		if (!(displayMin == 0.0 && displayMax == 0.0)) {
			int type = imp.getType();
			ImageProcessor ip = imp.getProcessor();
			if (type == ImagePlus.GRAY8 || type == ImagePlus.COLOR_256)
				ip.setMinAndMax(displayMin, displayMax);
			else if (type == ImagePlus.GRAY16 || type == ImagePlus.GRAY32) {
				if (ip.getMin() != displayMin || ip.getMax() != displayMax)
					ip.setMinAndMax(displayMin, displayMax);
			}
		}

		int stackSize = imp.getStackSize();
		if (stackSize > 1) {
			int channels = (int) getDouble(props, "channels");
			int slices = (int) getDouble(props, "slices");
			int frames = (int) getDouble(props, "frames");
			if (channels == 0)
				channels = 1;
			if (slices == 0)
				slices = 1;
			if (frames == 0)
				frames = 1;
			// IJ.log("setCalibration: "+channels+"  "+slices+"  "+frames);
			if (channels * slices * frames == stackSize) {
				imp.setDimensions(channels, slices, frames);
				if (getBoolean(props, "hyperstack"))
					imp.setOpenAsHyperStack(true);
			}
		}
	}

	private String getTltPath(String path) {
		return path.replace(".mrc", ".tlt");
	}

	/************** Helper methods to convert numbers extracted from Properties ****/
	private Double getNumber(Properties props, String key) {
		String s = props.getProperty(key);
		if (s != null) {
			try {
				return Double.valueOf(s);
			} catch (NumberFormatException e) {
			}
		}
		return null;
	}

	private double getDouble(Properties props, String key) {
		Double n = getNumber(props, key);
		return n != null ? n.doubleValue() : 0.0;
	}

	private boolean getBoolean(Properties props, String key) {
		String s = props.getProperty(key);
		return s != null && s.equals("true") ? true : false;
	}

	private void setStackDisplayRange(ImagePlus imp) {
		ImageStack stack = imp.getStack();
		double min = Double.MAX_VALUE;
		double max = -Double.MAX_VALUE;
		int n = stack.getSize();
		for (int i = 1; i <= n; i++) {

			ImageProcessor ip = stack.getProcessor(i);
			ip.resetMinAndMax();
			if (ip.getMin() < min)
				min = ip.getMin();
			if (ip.getMax() > max)
				max = ip.getMax();
		}
		imp.getProcessor().setMinAndMax(min, max);
		imp.updateAndDraw();
	}

	/**
	 * @return the resize
	 */
	public boolean isResize() {
		return resize;
	}

	/**
	 * @param resize
	 *            the resize to set
	 */
	public void setResize(boolean resize) {
		this.resize = resize;
	}
}
