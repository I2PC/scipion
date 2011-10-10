package browser.imageitems;

import browser.imageitems.tableitems.AbstractGalleryImageItem;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileInfo;
import ij.process.FloatProcessor;
import ij.process.StackStatistics;
import ij.util.Tools;
import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;
import xmipp.ImageDouble;
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
public class ImageConverter {

    public static ImagePlus convertToImageJ(ImageDouble image, String title) {
        return convertToImageJ(image, title, true);
    }

    public static ImagePlus convertToImageJ(ImageDouble image, String title, boolean useLogarithm) {
        if (image.isPSD()) {
            image.convertPSD(useLogarithm);
        }

        int w = image.getXsize();
        int h = image.getYsize();
        int d = image.getZsize();
        long n = image.getNsize();

        ImageStack is = new ImageStack(w, h);

        for (long i = 0; i < n; i++) {
            for (int j = 0; j < d; j++) {
                double slice[] = image.getData(i, j);

                FloatProcessor processor = new FloatProcessor(w, h, slice);
                is.addSlice(String.valueOf(i), processor);
            }
        }

        ImagePlus imp = new ImagePlus(title, is);

        // Sets associated file info.
        File f = new File(title);
        FileInfo fi = new FileInfo();
        fi.directory = f.getParent();
        fi.fileName = f.getName();
        imp.setFileInfo(fi);

        return imp;
    }

    /*    public static ImagePlus convertToImageJ(Projection projection, String title) {
    int w = projection.getXsize();
    int h = projection.getYsize();
    
    FloatProcessor processor = new FloatProcessor(w, h, projection.getData());
    
    return new ImagePlus(title, processor);
    }*/
    public static ImagePlus convertToImageJ(MetaData md) {
        LinkedList<String> missing = new LinkedList<String>();
        ImagePlus imp = null;

        if (md.containsLabel(MDLabel.MDL_IMAGE)) {
            ImageStack is = null;

            long ids[] = md.findObjects();

            for (long id : ids) {
                String filename = md.getValueString(MDLabel.MDL_IMAGE, id, true);

                try {
                    ImageDouble image = new ImageDouble(filename);
                    ImagePlus slice = ImageConverter.convertToImageJ(image, filename);

                    if (is == null) {
                        is = new ImageStack(slice.getWidth(), slice.getHeight());
                    }

                    is.addSlice(filename, slice.getProcessor());
                } catch (Exception ex) {
                    ex.printStackTrace(System.err);
                    missing.add(filename);
                }
            }

            imp = new ImagePlus(md.getFilename(), is);

            // Sets associated file info.
            File f = new File(md.getFilename());
            FileInfo fi = new FileInfo();
            fi.directory = f.getParent();
            fi.fileName = f.getName();
            imp.setFileInfo(fi);

            // Normalize by default
            normalizeStack(imp);
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

    public static ImagePlus convertToImageJ(ArrayList<AbstractGalleryImageItem> items) {
        ImageStack is = null;

        for (int i = 0; i < items.size(); i++) {
            AbstractGalleryImageItem item = items.get(i);

            if (item.isEnabled() && item.exists()) {
                ImagePlus ipslice = item.getImagePlus();

                if (is == null) {
                    is = new ImageStack(ipslice.getWidth(), ipslice.getHeight());
                }

                is.addSlice(ipslice.getTitle(), ipslice.getProcessor());
            }
        }

        return new ImagePlus("", is);
    }

    public static void normalizeStack(ImagePlus imp) {
        StackStatistics ss = new StackStatistics(imp);
        imp.getProcessor().setMinAndMax(ss.min, ss.max);
    }

    public static void revert(ImagePlus imp, String path) throws Exception {
        ImageDouble image = new ImageDouble(path);
        int w = image.getXsize();
        int h = image.getYsize();
        int d = image.getZsize();
        long n = image.getNsize();

        int sliceSize = w * h;
        int imageSize = sliceSize * d;
        double data[] = image.getData();
        double slice[] = new double[sliceSize];
        ImageStack is = new ImageStack(w, h);

        for (int i = 0; i < n; i++) {
            int offset = i * imageSize;
            for (int j = 0; j < d; j++) {
                System.arraycopy(data, offset + j * sliceSize, slice, 0, sliceSize);

                FloatProcessor processor = new FloatProcessor(w, h, slice);
                is.addSlice(String.valueOf(i), processor);
            }
        }

        imp.setStack(is);
    }

    public static ImageDouble convertToXmipp(ImagePlus imp) {
        ImageDouble image = new ImageDouble();

        int w = imp.getWidth();
        int h = imp.getHeight();
        int d = imp.getStackSize();

        double data[] = new double[w * h * d];
        for (int i = 0; i < d; i++) {
            float slice[] = (float[]) imp.getStack().getProcessor(i + 1).getPixels();
            System.arraycopy(Tools.toDouble(slice), 0, data, i * w * h, w * h);
        }
        try {
            image.setData(w, h, d, data);
        } catch (Exception ex) {
            ex.printStackTrace(System.err);
            IJ.error(ex.getMessage());
        }

        return image;
    }

    public static boolean saveImage(ImagePlus imp, String filename) {
        try {
            ImageDouble image = convertToXmipp(imp);

            image.write(filename);
            return true;
        } catch (Exception ex) {
            ex.printStackTrace(System.err);
            IJ.error(ex.getMessage());
        }

        return false;
    }
}
