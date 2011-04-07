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
import java.util.Enumeration;
import java.util.Iterator;

import xmipp.*;

public class TiltSeriesIO {
	// this class uses JNI mechanisms to call C++ code through the xmipp package
	// if an image is bigger than this threshold, resize it to this size
	// This optimization makes sense since for interactive use no more detail is needed
	public static Dimension resizeThreshold = new Dimension(400,400);

	public static boolean isImage(String path){
		// right now, identify file type by its extension
		// should check also uppercase like .SEL
		return path.endsWith(".mrc") || path.endsWith(".mrcs") || path.endsWith(".stk") || path.endsWith(".vol") || path.endsWith(".spi");
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
	
	public static void read(TomoData model) throws java.io.IOException,
			InterruptedException {

		String absolutePath = model.getFilePath();

		if ((absolutePath == null) || (absolutePath.equals("")))
			throw new IOException("Empty path");
		
		// initial assumption: the paths inside the selfile  need no traslation / rebuild
		boolean buildPath=false;
		
		model.readMetadata(absolutePath);
		// Xmipp_Tomo.debug(fn.getBaseName());
		
		long ids[]=model.getImageIds();
		for(long id:ids) {
			String fileName=model.getFilename(id);
			if(buildPath)
				fileName=buildAbsolutePath(absolutePath, fileName);
			ImageDouble image = null;
			if(fileName != null){
				try {
					image= readSlice(fileName, null, false);
				} catch (IOException ex) {
					// maybe the exception was due to the need of absolute paths
					// give a second try with absolute path
					buildPath = true;
					fileName=buildAbsolutePath(absolutePath, fileName);
					// another option would be to change the current working directory, but Java does not allow for it
					image = readSlice(fileName, null, false);
					
				}
				if(image != null)
					postReadSlice(model, image);
			}
		}

		postReadStack(model);
		
		if(! isSelFile(absolutePath)){
			// read tilt angles
			String tltFilePath = getTiltFilePath(absolutePath);
			model.readMetadata(tltFilePath);
		}
	}

		// TODO: bug - write fails (the written file cannot be read back again)
	// TODO: write details (below)
	// First show dialog to choose file type to save (primarily sel vs stack)
	// Sel: ideally ask for sel and stack file paths (defaults to same paths, obviously warning about overwriting in both cases)
	//      right now use same base and change extension only - for instance, pp.sel & pp.stk
	//      adjust stack filenames in metadata if necessary
	// Stack: ask only for stack file path (same overwrite warning), then save both stack and tlt
	public static void write(TomoData model) {
		String path = model.getFilePath();
		Xmipp_Tomo.debug("writing " + path);
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
			Xmipp_Tomo.debug("write", ex);
		}
	}
	
	private static void writeStack(String absolutePath, TomoData model) throws Exception{	
		ImageDouble img = convert(model.getImage());
		img.setFilename(absolutePath);
		img.write(absolutePath);
	}

	// TODO: - CURRENT - writeSel
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

	private static ImageDouble readSlice(String path, String base, boolean headerOnly)
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
	 * 
	 * @param model where the postprocessed data is saved
	 * @param img
	 * @depends readImage
	 */
	private static void postReadSlice(TomoData model, ImageDouble img) {
		ImagePlus image = convert(img);

		model.setOriginalHeight(img.getYsize());
		model.setOriginalWidth(img.getXsize());
		model.updateMinMax(image);

		// resize/scale - use an aux image processor for all projections
		ImageProcessor ipresized = null, ip = image.getProcessor();

		if (shouldResize(img.getXsize(),img.getYsize())) {
			ipresized = ip.resize(resizeThreshold.width,resizeThreshold.height);
			model.addProjection(ipresized);
			model.setResized(true);
		} else
			model.addProjection(ip);
	}
	
	/**
	 * Actions required after reading all the images of a stack
	 * 
	 * @param model
	 */
	private static void postReadStack(TomoData model) {
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
			Xmipp_Tomo.debug("Null image");
			return null;
		}
		// Creates image
		int width=img.getWidth(), height=img.getHeight(), numberOfProjections=img.getStackSize();
		ImageDouble imageDouble=new ImageDouble();

		int imageSize=width*height;
		double data[]=new double[imageSize*numberOfProjections];
		int i=0;
		for(int p=1;p<=numberOfProjections;p++){
			FloatProcessor projection=(FloatProcessor) img.getStack().getProcessor(p);
			float [] pixels =(float[]) projection.getPixels();
			for (int j = 0; j < imageSize; j++)
					data[i++]=pixels[j]; 
		}
		
		try{
			imageDouble.setData(width, height, numberOfProjections, data);
		}catch (Exception ex){
			Xmipp_Tomo.debug("convert ImagePlus->ImageDouble", ex);
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
			Xmipp_Tomo.debug("readTiltAngles", ex);
		}
	}

	/**
	 * TODO: writeTiltAngles - test
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
	
	public static void test1(String filepath){
		try {
			for(int i=1;i<=3; i++){
				ImageDouble image = readSlice(filepath, null, false, i);
				String firstFilePath=changeExtension(filepath, "" + i + ".xmp");
				image.write(firstFilePath);
			}
		} catch (Exception ex){
			Xmipp_Tomo.debug("test1",ex);
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
			Xmipp_Tomo.debug("test1",ex);
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
