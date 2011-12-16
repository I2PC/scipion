package xmippij;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileInfo;
import ij.process.ImageProcessor;
import ij.process.ByteProcessor;
import ij.process.ShortProcessor;
import ij.process.FloatProcessor;
import ij.process.StackStatistics;
import java.io.File;
import java.util.LinkedList;
import xmipp.Filename;
import xmipp.ImageGeneric;
import xmipp.MDLabel;
import xmipp.MetaData;


/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippImageConverter {

    public static ImagePlus loadImage(String filename) throws Exception {
        return loadImage(filename, 100);
    }

    public static ImagePlus loadImage(String filename, int zoom) throws Exception {
        ImageGeneric image = new ImageGeneric(filename);

        double scale = zoom / 100.0;
        int w = (int) (image.getXDim() * scale);
        int h = (int) (image.getYDim() * scale);

        return convertToImageJ(image, w, h);
    }

    public static ImagePlus convertToImageJ(ImageGeneric image) throws Exception {
        return convertToImageJ(image, ImageGeneric.ALL_SLICES);
    }

    public static ImagePlus convertToImageJ(ImageGeneric image, int nslice) throws Exception {
        return convertToImageJ(image, image.getXDim(), image.getYDim(), nslice);
    }

    public static ImagePlus convertToImageJ(ImageGeneric image, long nimage) throws Exception {
        return convertToImageJ(image, image.getXDim(), image.getYDim(), nimage);
    }

    public static ImagePlus convertToImageJ(ImageGeneric image, int width, int height) throws Exception {
        return convertToImageJ(image, width, height, ImageGeneric.ALL_SLICES, ImageGeneric.ALL_IMAGES);
    }

    public static ImagePlus convertToImageJ(ImageGeneric image, int width, int height, int nslice) throws Exception {
        return convertToImageJ(image, width, height, nslice, ImageGeneric.ALL_IMAGES);
    }

    public static ImagePlus convertToImageJ(ImageGeneric image, int width, int height, long nimage) throws Exception {
        return convertToImageJ(image, width, height, ImageGeneric.ALL_SLICES, nimage);
    }

    public static ImagePlus convertToImageJ(ImageGeneric image, int width, int height, int nslice, long nimage) throws Exception {
        ImageStack is = new ImageStack(width, height);

        ProcessorCreator pc = createProcessorCreator(image);

        boolean retrieveAllImages = nimage == ImageGeneric.ALL_IMAGES;
        boolean retrieveAllSlices = nslice == ImageGeneric.ALL_SLICES;
        long n = retrieveAllImages ? ImageGeneric.FIRST_IMAGE : nimage;
        long N = image.getNDim();

        for (; n <= N; n++) {
            image.read(width, height, n);

            if (image.isPSD()) {
                image.convertPSD(image.getUseLogarithm());
            }

            // Read volume.
            int slice = retrieveAllSlices ? ImageGeneric.FIRST_SLICE : nslice;
            for (; slice <= image.getZDim(); slice++) {
                ImageProcessor processor = pc.getProcessor(image, slice);
                ImagePlus imp = new ImagePlus("", processor);
                normalizeImagePlus(imp);
                is.addSlice("", imp.getProcessor());

                // If just one image, breaks loop.
                if (!retrieveAllSlices) {
                    break;
                }
            }

            // If just one image, breaks loop.
            if (!retrieveAllImages) {
                break;
            }
        }

        return buildImagePlus(image.getFilename(), is);
    }

    /**
     * Converts an ImageGeneric to ImageJ
     * WARNING!: Use this when converting images from memory.
     * @param image
     * @return
     * @throws Exception 
     */
    public static ImagePlus convertImageGenericToImageJ(ImageGeneric image) throws Exception {
        return convertImageGenericToImageJ(image, ImageGeneric.ALL_SLICES);
    }

    /**
     * Converts an ImageGeneric to ImageJ
     * WARNING!: Use this when converting images from memory.
     * @param image
     * @param nslice
     * @return
     * @throws Exception 
     */
    static ImagePlus convertImageGenericToImageJ(ImageGeneric image, int nslice) throws Exception {
        int width = image.getXDim();
        int height = image.getYDim();
        int depth = image.getZDim();
        ImageStack is = new ImageStack(width, height);

        ProcessorCreator pc = createProcessorCreator(image);

        boolean retrieveAllSlices = nslice == ImageGeneric.ALL_SLICES;

        // Read volume.      
        int slice = retrieveAllSlices ? ImageGeneric.FIRST_SLICE : nslice;
        for (; slice <= depth; slice++) {
            ImageProcessor processor = pc.getProcessor(image, slice);
            is.addSlice("", processor);

            // If just one image, breaks loop.
            if (!retrieveAllSlices) {
                break;
            }
        }

        ImagePlus imp = buildImagePlus(image.getFilename(), is);
        normalizeImagePlus(imp);

        return imp;
    }

    public static ImagePlus convertToImageJ(MetaData md) throws Exception {
        LinkedList<String> missing = new LinkedList<String>();
        ImagePlus imp = null;

        if (md.containsLabel(MDLabel.MDL_IMAGE)) {
            ImageStack is = null;

            long ids[] = md.findObjects();

            for (long id : ids) {
                String filename = md.getValueString(MDLabel.MDL_IMAGE, id, true);

                try {
                    String name = Filename.getFilename(filename);
                    long n = Filename.getNimage(filename);

                    ImageGeneric image = new ImageGeneric(filename);
                    ImagePlus slice = convertToImageJ(image, n);

                    if (is == null) {
                        is = new ImageStack(slice.getWidth(), slice.getHeight());
                    }

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
    }

    public static ImagePlus buildImagePlus(String filename, ImageStack is) {
        ImagePlus imp = new ImagePlus(filename, is);
        // Sets associated file info.
        File f = new File(filename);
        FileInfo fi = new FileInfo();
        String absPath = f.getAbsolutePath();
        fi.directory = absPath.substring(0, absPath.lastIndexOf(File.separator));
        fi.fileName = f.getName();

        imp.setFileInfo(fi);

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

    public static ImageGeneric convertToXmipp(ImagePlus imp) throws Exception {
        ImageGeneric image = new ImageGeneric();

        DataSetter ds = createDataSetter(image, imp.getType());

        image.resize(imp.getWidth(), imp.getHeight(), imp.getStackSize());

        ImageStack is = imp.getStack();
        for (int nslice = ImageGeneric.FIRST_SLICE; nslice <= is.getSize(); nslice++) {
            ds.setArray(image, is.getProcessor(nslice).getPixels(), nslice);
        }

        return image;
    }

    public static boolean saveImage(ImagePlus imp, String filename) throws Exception {
        ImageGeneric image = new ImageGeneric();

        int width = imp.getWidth();
        int height = imp.getHeight();
        int depth = imp.getStackSize();
        int image_index = 1;

        boolean storeAllImages = Filename.hasStackExtension(filename);
        boolean storeAllSlices = !storeAllImages;// = Filename.hasVolumeExtension(filename);

        // @TODO Add support for hyperstacks and remove this "if".
        if (storeAllImages) {
            image_index = depth;
            depth = 1;
        }

        DataSetter ds = createDataSetter(image, imp.getType());

        for (long nimage = ImageGeneric.FIRST_IMAGE; nimage <= image_index; nimage++) {
            image.mapFile2Write(width, height, depth, filename, nimage);

            // Store volume.
            for (int nslice = ImageGeneric.FIRST_SLICE; nslice <= depth; nslice++) {
                Object data = imp.getStack().getProcessor((int) nimage * nslice).getPixels();

                // image.setData...
                ds.setArray(image, data, nslice);

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

    static ProcessorCreator createProcessorCreator(ImageGeneric image) throws Exception {
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

    static DataSetter createDataSetter(ImageGeneric image, int dataType) throws Exception {

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

    public abstract ImageProcessor getProcessor(ImageGeneric image, int slice) throws Exception;
}

class ProcessorCreatorByte extends ProcessorCreator {

    @Override
    public ImageProcessor getProcessor(ImageGeneric image, int slice) throws Exception {
        return new ByteProcessor(image.getXDim(), image.getYDim(), image.getArrayByte(slice), null);
    }
}

class ProcessorCreatorShort extends ProcessorCreator {

    @Override
    public ImageProcessor getProcessor(ImageGeneric image, int slice) throws Exception {
        return new ShortProcessor(image.getXDim(), image.getYDim(), image.getArrayShort(slice), null);
    }
}

class ProcessorCreatorFloat extends ProcessorCreator {

    @Override
    public ImageProcessor getProcessor(ImageGeneric image, int slice) throws Exception {
        return new FloatProcessor(image.getXDim(), image.getYDim(), image.getArrayFloat(slice), null);
    }
}

abstract class DataSetter {

    public abstract void setArray(ImageGeneric image, Object data, int slice) throws Exception;
}

class DataSetterByte extends DataSetter {

    @Override
    public void setArray(ImageGeneric image, Object data, int slice) throws Exception {
        image.setArrayByte((byte[]) data, slice);
    }
}

class DataSetterShort extends DataSetter {

    @Override
    public void setArray(ImageGeneric image, Object data, int slice) throws Exception {
        image.setArrayShort((short[]) data, slice);
    }
}

class DataSetterFloat extends DataSetter {

    @Override
    public void setArray(ImageGeneric image, Object data, int slice) throws Exception {
        image.setArrayFloat((float[]) data, slice);
    }
}
