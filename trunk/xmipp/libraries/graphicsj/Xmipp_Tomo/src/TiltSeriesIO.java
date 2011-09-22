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
 * @deprecated - use ImageDouble and Metadata directly
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
import java.util.Enumeration;


import xmipp.*;

/**
 * @deprecated
 * TODO: remove when not needed anymore...
 * @author jcuenca
 *
 */
public class TiltSeriesIO {
	// this class uses JNI mechanisms to call C++ code through the xmipp package

	/**
	 * @deprecated
	 */
	public static Dimension resizeThreshold = new Dimension(400,400);

	public static boolean isImage(String path){
		// right now, identify file type by its extension
		// should check also uppercase like .SEL
		return path.endsWith(".mrc") || path.endsWith(".mrcs") || path.endsWith(".stk") || path.endsWith(".vol") || path.endsWith(".spi");
	}
	
	public static boolean isAbsolute(String path){
		int indexAt = path.indexOf("@");
		if(indexAt >= 0)
			return path.charAt(indexAt+1) == '/';
		else{
			File test=new File(path);
			return test.isAbsolute();
		}
	}
	
	public static boolean isSelFile(String path){
		// right now, identify file type by its extension
		// should check also uppercase like .SEL
		return path.endsWith(".sel");
	}
	
	public static boolean isTltFile(String path){
		// right now, identify file type by its extension
		// should check also uppercase like .SEL
		return path.endsWith(".tlt");
	}
	
	// fixed: fails with mrcs due to memory leak (bad alloc: coreAllocateReuse tries to alloc X*Y*Z*N bytes,
	// with N being the number of slices, even though we are requesting to read only 1)
	public static void read(TomoData model) throws java.io.IOException,
			InterruptedException {

		String absolutePath = model.getFilePath();

		if ((absolutePath == null) || (absolutePath.equals("")))
			throw new IOException("Empty path");
		
		// initial assumption: the paths inside the selfile  need no traslation / rebuild
		boolean buildPath=false;
		
		// Xmipp_Tomo.debug("Path: "+ absolutePath + ". Build: " + buildPath);
		model.readMetadata(absolutePath);
		
		long ids[]=model.getStackIds();
		for(long id:ids) {
			String fileName=model.getFilePath(id);
			if(buildPath)
				fileName=buildAbsolutePath(absolutePath, fileName);
			ImageDouble image = null;
			if(fileName != null){
				try {
					image= readSlice(fileName, null, false);
				} catch (IOException ex) {
					// maybe the exception was due to the need of absolute paths
					// give a second try with absolute path
					Logger.debug("TiltseriesIO.read", ex);
					if(! isAbsolute(fileName)){
						buildPath = true;
						fileName=buildAbsolutePath(absolutePath, fileName);
						// another option would be to change the current working directory, but Java does not allow for it
						image = readSlice(fileName, null, false);
					}
				}
				if(image != null)
					postReadSlice(model, image);
			}
		}

		postReadStack(model);
		
	}

	// fixed: bug - write fails? spider and derivatives fixed, test with other formats (MRC / MRCS)
	public static void write(TomoData model) {
		String path = model.getFilePath();
		Logger.debug("writing " + path);
		try{
			if(isImage(path)){
				writeStack(path,model);
				String tltFilePath = getTiltFilePath(path);
				writeTiltAngles(tltFilePath, model);
			}else if(isSelFile(path)){
				String stackPath=changeExtension(path, ".stk");
				writeSel(path,stackPath,model);
			}
				
		}catch (Exception ex){
			Logger.debug("write", ex);
		}
	}
	
	/**
	 * @deprecated
	 * @param absolutePath
	 * @param model
	 * @throws Exception
	 */
	private static void writeStack(String absolutePath, TomoData model) throws Exception{	
		ImageDouble img = convert(model.getImage());
		img.setFilename(absolutePath);
		// @see rwSpider.cpp -- WRITE_OVERWRITE required for maxim to be updated (readSpider gets nDim from maxim)
		img.write(absolutePath, 0, true, ImageWriteMode.WRITE_OVERWRITE, CastWriteMode.CW_CAST);
	}

	
	private static void writeSel(String selAbsolutePath, String stackAbsolutePath, TomoData model) throws Exception{

		model.getMetadata().write(selAbsolutePath);
	}


	

private static String buildAbsolutePath(String selFilePath, String path){
	// syntax of absolute paths: slice@/absolute/path/file.ext
	String baseDir=new File(selFilePath).getParent();
	int separatorIndex=path.indexOf("@");
	String slice = path.substring(0,separatorIndex+1);
	String fileName = path.substring(separatorIndex+1,path.length());
	return slice + baseDir + "/" + fileName;
}

	public static ImageDouble readSlice(String path, String base, boolean headerOnly)
			throws IOException {
		// watch out: slice = -1 has special meaning for image_base.read
		return readSlice(path, base, headerOnly, -10);
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
	 *            - slice number (1..N) (depends on Metadata convention for numbering slices)
	 * @return image as ImageDouble
	 * @throws IOException
	 */
	private static ImageDouble readSlice(String path, String base, boolean headerOnly,
			int slice) throws IOException {
		String filepath = path;
		if (base != null)
			filepath = base + "/" + path;

		// Xmipp_Tomo.debug(filepath);

		ImageDouble image = new ImageDouble();
		try{
			if(headerOnly)
				image.readHeader(filepath);
			else
				image.read(filepath, slice);
		}catch (Exception ex){
			throw new IOException(ex);
		}
		return image;
	}

	/**
	 * Actions required after reading 1 image
	 * @deprecated
	 * @param model where the postprocessed data is saved
	 * @param img
	 * @depends readImage
	 */
	public static void postReadSlice(TomoData model, ImageDouble img) {
		ImagePlus image = convert(img);

		model.setOriginalHeight(img.getYsize());
		model.setOriginalWidth(img.getXsize());
		model.updateMinMax(image);

		// resize/scale - use an aux image processor for all projections
		ImageProcessor ipresized = null, ip = image.getProcessor();

		if (shouldResize(img.getXsize(),img.getYsize())) {
			// ipresized = ip.resize(resizeThreshold.width,resizeThreshold.height);
			model.addProjection(ipresized);
			model.setResized(true);
		} else
			model.addProjection(ip);
	}
	
	/**
	 * Actions required after reading all the images of a stack
	 * @deprecated
	 * @param model
	 */
	public static void postReadStack(TomoData model) {
		model.lastImageLoaded();
	}

	
	private static ImagePlus convert(ImageDouble img) {
		double [] imageData = img.getData();
		if (imageData == null)
			return null;
		
		// Create image
		int width=img.getXsize(), height=img.getYsize();
		FloatProcessor convertedImageProcessor = new FloatProcessor(width,height);
		int i=0;
		for (int y = 0; y < height; y++)
			for (int x = 0; x < width; x++)	
					convertedImageProcessor.setf(x, y, (float) imageData[i++]);

		// normalize the image - done by default when calling FloatProcessor(array), here is our responsibility
		convertedImageProcessor.findMinAndMax();
		
		ImagePlus imagePlus = new ImagePlus("ImageDouble", convertedImageProcessor);
		return imagePlus;
	}
	
	private static ImageDouble convert(ImagePlus img) {
		if(img == null){
			Logger.debug("Null image");
			return null;
		}
		// Creates image 
		int width=img.getWidth(), height=img.getHeight(), numberOfProjections=img.getStackSize();
		ImageDouble imageDouble=new ImageDouble();

		int imageSize=width*height;
		double data[]=new double[imageSize*numberOfProjections];
		int i=0;
		// ImagePlus does not return the whole stack as a 3D array. 
		// Maybe this loop will be faster if block operations are used - like replacing the inner loop with a memcpy
		// Or, storing the float pixels array directly into the ImageDouble, with a new method like setProjection(float [])
		for(int p=1;p<=numberOfProjections;p++){
			FloatProcessor projection=(FloatProcessor) img.getStack().getProcessor(p);
			float [] pixels =(float[]) projection.getPixels();
			for (int j = 0; j < imageSize; j++)
					data[i++]=pixels[j]; 
		}
		
		try{
			// for tomograms we need N instead of Z, so Z = 1
			imageDouble.setData(width, height, 1, numberOfProjections, data);
		}catch (Exception ex){
			Logger.debug("convert ImagePlus->ImageDouble", ex);
		}
		return imageDouble;
	}

	// TODO: Metadata can also import columns from text files
	public static void readTiltAngles(String tltFilePath, MetaData md) {
		// .tlt syntax: one float angle per line, stored as text
		try{		
			BufferedReader brin = new BufferedReader(new FileReader(tltFilePath));
	
			String line = null;
			long i=1;
			while ((line = brin.readLine()) != null) {
				md.setValueDouble(MDLabel.MDL_ANGLETILT, new Double(line), i++);
			}
		}catch (IOException ex){
			// tilt file missing - do nothing
			// Xmipp_Tomo.debug("readTiltAngles", ex);
		}
	}

	/**
	 * TODO: writeTiltAngles - test to ensure the .tlt file is written correctly
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
			Enumeration<Double> tiltAngles = model.getTiltAngles();
			while (tiltAngles.hasMoreElements()){
				output.println(tiltAngles.nextElement());
			}

		} catch (FileNotFoundException ex) {
			Logger.debug("writeTlt - file not found");
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

/**
 * @deprecated
 * @param width
 * @param height
 * @return
 */
	public static boolean shouldResize(int width, int height){
			return (width > resizeThreshold.width ) || (height > resizeThreshold.height);
	}
	
	public static void test1(String filepath){
		try {
			for(int i=1;i<=3; i++){
				ImageDouble image = readSlice(filepath, null, false, i);
				String firstFilePath=changeExtension(filepath, "" + i + ".xmp");
				image.write(firstFilePath);
			}
		} catch (Exception ex){
			Logger.debug("test1",ex);
		}
	}
	
	public static void test2(String filepath){
		try {
			TomoData model=new TomoData(filepath);
			TiltSeriesIO.read(model);
			String extension=filepath.substring(filepath.lastIndexOf("."),filepath.length());
			String copyFilePath=changeExtension(filepath, "copy"+extension);
			model.setFile(copyFilePath);
			TiltSeriesIO.write(model);
		} catch (Exception ex){
			Logger.debug("test1",ex);
		}
	}

	public static void main(String[] args){
		if(args.length < 1){
			System.out.println("Please specify the data file for tests");
			System.exit(1);
		}
		String filepath=args[0];
		test1(filepath);
	}

}
