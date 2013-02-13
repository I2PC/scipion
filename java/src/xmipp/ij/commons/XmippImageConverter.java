
/***************************************************************************
 * Authors:     Juanjo Vega
 * 		 		J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
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

package xmipp.ij.commons;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileInfo;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import ij.process.StackStatistics;

import java.io.File;
import java.util.LinkedList;

import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MetaData;

public class XmippImageConverter {

	public static ImagePlus loadImage(String filename) throws Exception {
		return loadImage(filename, 100);
	}

	public static ImagePlus loadImage(String filename, int zoom)
			throws Exception {
		ImageGeneric image = new ImageGeneric(filename);

		double scale = zoom / 100.0;
		int w = (int) (image.getXDim() * scale);
		int h = (int) (image.getYDim() * scale);

		ImagePlus imp;
		/*
		 * if (Filename.hasPrefix(filename)) { String name =
		 * Filename.getFilename(filename); long nimage =
		 * Filename.getNimage(filename);
		 * 
		 * image = new ImageGeneric(name); imp = convertToImageJ(image, nimage);
		 * } else {
		 */
		long select_image = Filename.getNimage(filename);
		imp = readToImagePlus(image, w, h, select_image);
		// }

		return imp;
	}

	/**
	 * These set of functions will read from disk an ImageGeneric that had read
	 * the header previously. The result is an ImagePlus
	 */
	public static ImagePlus readToImagePlus(ImageGeneric image)
			throws Exception {
		return readToImagePlus(image, ImageGeneric.ALL_SLICES);
	}

	public static ImagePlus readToImagePlus(ImageGeneric image, int nslice)
			throws Exception {
		return readToImagePlus(image, image.getXDim(), image.getYDim(), nslice);
	}

	public static ImagePlus readToImagePlus(ImageGeneric image, long nimage)
			throws Exception {
		return readToImagePlus(image, image.getXDim(), image.getYDim(), nimage);
	}

	public static ImagePlus readToImagePlus(ImageGeneric image, int width,
			int height) throws Exception {
		return readToImagePlus(image, width, height, ImageGeneric.ALL_SLICES,
				ImageGeneric.ALL_IMAGES);
	}

	public static ImagePlus readToImagePlus(ImageGeneric image, int width,
			int height, int nslice) throws Exception {
		return readToImagePlus(image, width, height, nslice,
				ImageGeneric.ALL_IMAGES);
	}

	public static ImagePlus readToImagePlus(ImageGeneric image, int width,
			int height, long nimage) throws Exception {
		return readToImagePlus(image, width, height, ImageGeneric.ALL_SLICES,
				nimage);
	}

	/**
	 * In this function is suposed the the ImageGeneric is read only the header
	 * and the data will be read from disk
	 */
	public static ImagePlus readToImagePlus(ImageGeneric image, int width,
			int height, int slice, long select_image) throws Exception {
		ImageStack is = new ImageStack(width, height);
		ProcessorCreator pc = createProcessorCreator(image);
		long lastImage = select_image;

		if (select_image == ImageGeneric.ALL_IMAGES) {
			select_image = ImageGeneric.FIRST_IMAGE;
			lastImage = image.getNDim();
		}
		for (; select_image <= lastImage; select_image++) {
			image.read(width, height, select_image);
			addSlicesToStack(image, ImageGeneric.FIRST_IMAGE, slice, is, pc);
		}

		return buildImagePlus(image.getFilename(), is);
	}

	/**
	 * Converts an ImageGeneric to ImageJ WARNING!: Use this when converting
	 * images from memory.
	 */
	public static ImagePlus convertToImagePlus(ImageGeneric image)
			throws Exception {
		return convertToImagePlus(image, ImageGeneric.ALL_IMAGES,
				ImageGeneric.ALL_SLICES);
	}

	/**
	 * Converts an ImageGeneric to ImageJ WARNING!: Use this when converting
	 * images from memory.
	 */
	public static ImagePlus convertToImagePlus(ImageGeneric image,
			long select_image) throws Exception {
		return convertToImagePlus(image, select_image, ImageGeneric.ALL_SLICES);
	}

	public static ImagePlus convertToImagePlus(ImageGeneric image,
			long select_image, int slice) throws Exception {
		int width = image.getXDim();
		int height = image.getYDim();

		ImageStack is = new ImageStack(width, height);
		ProcessorCreator pc = createProcessorCreator(image);
		long lastImage = select_image;

		if (select_image == ImageGeneric.ALL_IMAGES) {
			select_image = ImageGeneric.FIRST_IMAGE;
			lastImage = image.getNDim();
		}
		for (; select_image <= lastImage; select_image++) {
			addSlicesToStack(image, select_image, slice, is, pc);
		}
		return buildImagePlus(image.getFilename(), is);
	}

	/** Internal function to add the slices to ImageStack to build the ImagePlus */
	private static void addSlicesToStack(ImageGeneric image, long select_image,
			int slice, ImageStack is, ProcessorCreator pc) throws Exception {
		if (image.isPSD()) {
			image.convertPSD(image.getUseLogarithm());
		}
		int lastSlice = slice;
		if (slice == ImageGeneric.ALL_SLICES) {
			slice = ImageGeneric.FIRST_SLICE;
			lastSlice = image.getZDim();
		}
		for (; slice <= lastSlice; slice++)
			is.addSlice("", pc.getProcessor(image, select_image, slice));

	}

	public static ImagePlus readMdRowToImagePlus(String fn, MetaData md,
			long objId, int width, int height, boolean useGeo, boolean wrap)
			throws Exception {
		ImageGeneric image = new ImageGeneric();
		if (useGeo)
			image.readApplyGeo(fn, md, objId, width, height, wrap);
		else
			image.read(fn, width, height);
		return convertToImagePlus(image);
	}

	public static ImagePlus readMdRowToImagePlus(String fn, MetaData md,
			long objId, boolean useGeo, boolean wrap) throws Exception {
		ImageGeneric image = new ImageGeneric();
		if (useGeo)
			image.readApplyGeo(fn, md, objId, wrap);
		else
			image.read(fn);
		return convertToImagePlus(image);
	}

	/** Read the entries on the metadata and create an ImagePlus */
	public static ImagePlus readMetadataToImagePlus(int label, MetaData md,
			boolean useGeo, boolean wrap) throws Exception {
		LinkedList<String> missing = new LinkedList<String>();
		ImagePlus imp = null;

		if (md.containsLabel(label)) {
			ImageStack is = null;
			long ids[] = md.findObjects();

			for (long id : ids) {
				String filename = md.getValueString(label, id, true);

				try {
					ImagePlus slice = readMdRowToImagePlus(filename, md, id, useGeo, wrap);

					if (is == null)
						is = new ImageStack(slice.getWidth(), slice.getHeight());
					is.addSlice(filename, slice.getProcessor());
				} catch (Exception ex) {
					ex.printStackTrace(System.err);
					missing.add(filename);
				}
			}
			imp = buildImagePlus(md.getFilename(), is);
		}
		// Tells user about missing files.
		if (!missing.isEmpty()) {
			String message = "There are missing files:\n";
			for (int i = 0; i < missing.size(); i++) {
				message += missing.get(i) + "\n";
			}
			IJ.error(message);
		}

		return imp;
	}//function readMetadataToImagePlus

	public static ImagePlus buildImagePlus(String filename, ImageStack is) {
		ImagePlus imp = new ImagePlus(filename, is);

		if (filename != null) {
			// Sets associated file info.
			File f = new File(filename);
			FileInfo fi = new FileInfo();
			String absPath = f.getAbsolutePath();
			fi.directory = absPath.substring(0,
					absPath.lastIndexOf(File.separator));
			fi.fileName = f.getName();

			imp.setFileInfo(fi);
		}

		return imp;
	}

	public static void normalizeImagePlus(ImagePlus imp) {
		StackStatistics ss = new StackStatistics(imp);
		imp.getProcessor().setMinAndMax(ss.min, ss.max);
	}

	public static void revert(ImagePlus imp, String path) throws Exception {
		ImagePlus imp2 = XmippImageConverter.loadImage(path);
		imp2.setTitle("Revert!!: " + System.currentTimeMillis());
		imp2.show();
		imp.setStack(imp.getTitle(), imp2.getImageStack());
	}

	/**
	 * This will convert an ImagePlus to ImageGeneric. NOTE: This can be used
	 * now with images or volumes
	 */
	public static ImageGeneric convertToImageGeneric(ImagePlus imp)
			throws Exception {
		ImageGeneric image = new ImageGeneric();
		DataSetter ds = createDataSetter(image, imp.getType());
		image.resize(imp.getWidth(), imp.getHeight(), imp.getStackSize());

		ImageStack is = imp.getStack();
		for (int nslice = ImageGeneric.FIRST_SLICE; nslice <= is.getSize(); nslice++) {
			ds.setArray(image, is.getProcessor(nslice).getPixels(),
					ImageGeneric.FIRST_IMAGE, nslice);
		}

		return image;
	}

	public static boolean writeImagePlus(ImagePlus imp, String filename)
			throws Exception {
		ImageGeneric image = new ImageGeneric();

		int width = imp.getWidth();
		int height = imp.getHeight();
		int depth = imp.getStackSize();
		int image_index = 1;

		boolean storeAllImages = Filename.hasStackExtension(filename);
		boolean storeAllSlices = !storeAllImages;// =
													// Filename.hasVolumeExtension(filename);

		// TODO: Add support for hyperstacks and remove this "if".
		if (storeAllImages) {
			image_index = depth;
			depth = 1;
		}

		DataSetter ds = createDataSetter(image, imp.getType());

		for (long nimage = ImageGeneric.FIRST_IMAGE; nimage <= image_index; nimage++) {
			image.mapFile2Write(width, height, depth, filename, nimage);

			// Store volume.
			for (int nslice = ImageGeneric.FIRST_SLICE; nslice <= depth; nslice++) {
				Object data = imp.getStack()
						.getProcessor((int) nimage * nslice).getPixels();

				// image.setData...
				ds.setArray(image, data, ImageGeneric.FIRST_IMAGE, nslice);

				// If just one image, breaks loop.
				if (!storeAllSlices) {
					break;
				}
			}

			image.write(filename);

			// If just one image, breaks loop.
			if (!storeAllImages) {
				break;
			}
		}

		return true;
	}

	static ProcessorCreator createProcessorCreator(ImageGeneric image)
			throws Exception {
		ProcessorCreator pc = null;
		switch (image.getDataType()) {
		case ImageGeneric.Float:
		case ImageGeneric.Double:
			pc = new ProcessorCreatorFloat();
			break;
		case ImageGeneric.Short:
		case ImageGeneric.UShort:
			pc = new ProcessorCreatorShort();
			break;
		case ImageGeneric.SChar:
		case ImageGeneric.UChar:
			pc = new ProcessorCreatorByte();
			break;
		default:
			pc = new ProcessorCreatorFloat();
		}

		return pc;
	}

	static DataSetter createDataSetter(ImageGeneric image, int dataType)
			throws Exception {

		// Creates a DataSetter according to dataType.
		DataSetter ds = null;
		boolean setDataType = image.getDataType() == ImageGeneric.Unknown_Type;

		switch (dataType) {
		case ImagePlus.GRAY32:
			if (setDataType) {
				image.setDataType(ImageGeneric.Float);
			}
			ds = new DataSetterFloat();
			break;
		case ImagePlus.GRAY16:
			if (setDataType) {
				image.setDataType(ImageGeneric.UShort);
			}
			ds = new DataSetterShort();
			break;
		case ImagePlus.GRAY8:
			if (setDataType) {
				image.setDataType(ImageGeneric.UChar);
			}
			ds = new DataSetterByte();
			break;
		default:
			throw new Exception("Format not supported: " + dataType);
		}

		return ds;
	}
}

abstract class ProcessorCreator {

	public abstract ImageProcessor getProcessor(ImageGeneric image,
			long select_image, int slice) throws Exception;
}

class ProcessorCreatorByte extends ProcessorCreator {

	@Override
	public ImageProcessor getProcessor(ImageGeneric image, long select_image,
			int slice) throws Exception {
		return new ByteProcessor(image.getXDim(), image.getYDim(),
				image.getArrayByte(select_image, slice), null);
	}
}

class ProcessorCreatorShort extends ProcessorCreator {

	@Override
	public ImageProcessor getProcessor(ImageGeneric image, long select_image,
			int slice) throws Exception {
		return new ShortProcessor(image.getXDim(), image.getYDim(),
				image.getArrayShort(select_image, slice), null);
	}
}

class ProcessorCreatorFloat extends ProcessorCreator {

	@Override
	public ImageProcessor getProcessor(ImageGeneric image, long select_image,
			int slice) throws Exception {
		return new FloatProcessor(image.getXDim(), image.getYDim(),
				image.getArrayFloat(select_image, slice), null);
	}
}

abstract class DataSetter {

	public abstract void setArray(ImageGeneric image, Object data,
			long select_image, int slice) throws Exception;
}

class DataSetterByte extends DataSetter {

	@Override
	public void setArray(ImageGeneric image, Object data, long select_image,
			int slice) throws Exception {
		image.setArrayByte((byte[]) data, select_image, slice);
	}
}

class DataSetterShort extends DataSetter {

	@Override
	public void setArray(ImageGeneric image, Object data, long select_image,
			int slice) throws Exception {
		image.setArrayShort((short[]) data, select_image, slice);
	}
}

class DataSetterFloat extends DataSetter {

	@Override
	public void setArray(ImageGeneric image, Object data, long select_image,
			int slice) throws Exception {
		image.setArrayFloat((float[]) data, select_image, slice);
	}
}
