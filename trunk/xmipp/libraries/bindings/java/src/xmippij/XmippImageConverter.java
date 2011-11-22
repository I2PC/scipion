package xmippij;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileInfo;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import ij.process.FloatProcessor;
import ij.process.StackStatistics;
import java.io.File;
import java.util.LinkedList;
import xmipp.ImageGeneric;
import xmipp.Projection;
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

    public static ImagePlus convertToImagej(String filename) throws Exception {
        ImageGeneric image = new ImageGeneric(filename);
        return convertToImagej(image, true);
    }

    public static ImagePlus convertToImagej(ImageGeneric image) throws Exception {
        return convertToImagej(image, true);
    }

    public static ImagePlus convertToImagej(ImageGeneric image, boolean useLogarithm) throws Exception {
//TODO: FIXME
        //        if (image.isPSD()) {
//            image.convertPSD(useLogarithm);
//        }

        String filename = image.filename;
        int w = image.xSize;
        int h = image.ySize;
        int d = image.zSize;
        long n = image.nSize;

        ImageStack is = new ImageStack(w, h);
        ProcessorCreator pc = null;
        switch (image.dataType) {
            case ImageGeneric.Float:
            case ImageGeneric.Double:
                pc = new FloatProcessorCreator();
                break;
            case ImageGeneric.SChar:
            case ImageGeneric.UChar:
                pc = new ByteProcessorCreator();
                break;
            case ImageGeneric.Short:
            case ImageGeneric.UShort:
                pc = new ShortProcessorCreator();
                break;
            default:
                pc = new FloatProcessorCreator();
            //throw new Exception("Unrecognized image datatype " + image.dataType);
        }

        for (long nimage = ImageGeneric.FIRST_IMAGE; nimage <= n; nimage++) {
            for (int nslice = ImageGeneric.FIRST_SLICE; nslice <= d; nslice++) {
                image.nSize = nimage;
                image.zSize = nslice;
                is.addSlice(String.valueOf((nimage * d) + nslice), pc.getProcessor(image));
            }
        }

        return getNormalizedImagePlus(filename, is);
    }

    public static ImagePlus convertToImagej(MetaData md)
            throws Exception {
        LinkedList<String> missing = new LinkedList<String>();
        ImagePlus imp = null;

        if (md.containsLabel(MDLabel.MDL_IMAGE)) {
            ImageStack is = null;

            long ids[] = md.findObjects();

            for (long id : ids) {
                String filename = md.getValueString(MDLabel.MDL_IMAGE, id, true);

                try {
                    ImagePlus slice = convertToImagej(filename);

                    if (is == null) {
                        is = new ImageStack(slice.getWidth(), slice.getHeight());
                    }

                    is.addSlice(filename, slice.getProcessor());
                } catch (Exception ex) {
                    ex.printStackTrace(System.err);
                    missing.add(filename);
                }
            }

            imp = getNormalizedImagePlus(md.getFilename(), is);
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

    public static ImagePlus convertToImagej(Projection projection, String title) {
        FloatProcessor processor = new FloatProcessor(projection.getXsize(), projection.getYsize(), projection.getData());
        return new ImagePlus(title, processor);
    }

    public static ImagePlus getNormalizedImagePlus(String filename, ImageStack is) {
        ImagePlus imp = new ImagePlus(filename, is);
        // Sets associated file info.
        File f = new File(filename);
        FileInfo fi = new FileInfo();
        fi.directory = f.getParent();
        fi.fileName = f.getName();
        imp.setFileInfo(fi);

        // Normalize by default
        StackStatistics ss = new StackStatistics(imp);
        imp.getProcessor().setMinAndMax(ss.min, ss.max);
        return imp;
    }

    public static void revert(ImagePlus imp, String path) throws Exception {
        ImageGeneric image = new ImageGeneric(path);
        ImagePlus imp2 = XmippImageConverter.convertToImagej(image);

        imp.setStack(imp.getTitle(), imp2.getImageStack());
//        int w = image.xSize;
//        int h = image.ySize;
//        int d = image.zSize;
//        long n = image.nSize;
//
//        int sliceSize = w * h;
//        int imageSize = sliceSize * d;
//        double data[] = image.getData();
//        double slice[] = new double[sliceSize];
//        ImageStack is = new ImageStack(w, h);
//
//        for (int i = 0; i < n; i++) {
//            int offset = i * imageSize;
//            for (int j = 0; j < d; j++) {
//                System.arraycopy(data, offset + j * sliceSize, slice, 0, sliceSize);
//
//                FloatProcessor processor = new FloatProcessor(w, h, slice);
//                is.addSlice(String.valueOf(i), processor);
//            }
//        }
//
//        imp.setStack(is);
    }

    public static ImageGeneric convertToXmipp(ImagePlus imp) {
        ImageGeneric image = new ImageGeneric();

//        int w = imp.getWidth();
//        int h = imp.getHeight();
//        int d = imp.getStackSize();
//
//        double data[] = new double[w * h * d];
//        for (int i = 0; i < d; i++) {
//            float slice[] = (float[]) imp.getStack().getProcessor(i + 1).getPixels();
//            System.arraycopy(Tools.toDouble(slice), 0, data, i * w * h, w * h);
//        }
//        try {
//            image.setData(w, h, d, data);
//        } catch (Exception ex) {
//            ex.printStackTrace(System.err);
//            IJ.error(ex.getMessage());
//        }

        return image;
    }

    public static boolean saveImage(ImagePlus imp, String filename) {
//        try {
//            ImageGeneric image = convertToXmipp(imp);
//
//            image.write(filename);
//            return true;
//        } catch (Exception ex) {
//            ex.printStackTrace(System.err);
//            IJ.error(ex.getMessage());
//        }

        return false;
    }
}

abstract class ProcessorCreator {

    public abstract ImageProcessor getProcessor(ImageGeneric image) throws Exception;
}

class FloatProcessorCreator extends ProcessorCreator {

    public ImageProcessor getProcessor(ImageGeneric image) throws Exception {
        return new FloatProcessor(image.xSize, image.ySize, image.getArrayFloat(), null);
    }
}

class ShortProcessorCreator extends ProcessorCreator {

    public ImageProcessor getProcessor(ImageGeneric image) throws Exception {
        return new ShortProcessor(image.xSize, image.ySize, image.getArrayShort(), null);
    }
}

class ByteProcessorCreator extends ProcessorCreator {

    public ImageProcessor getProcessor(ImageGeneric image) throws Exception {
        return new ij.process.ByteProcessor(image.xSize, image.ySize, image.getArrayByte(), null);
    }
}