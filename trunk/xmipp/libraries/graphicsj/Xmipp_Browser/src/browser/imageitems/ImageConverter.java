package browser.imageitems;

import browser.imageitems.tableitems.AbstractTableImageItem;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileInfo;
import ij.process.FloatProcessor;
import ij.process.StackStatistics;
import ij.util.Tools;
import java.io.File;
import java.util.LinkedList;
import java.util.Vector;
import xmipp.ImageDouble;
import xmipp.MDLabel;
import xmipp.MetaData;
import xmipp.Projection;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class ImageConverter {

    public static ImagePlus convertToImagej(ImageDouble image, String title) {
        return convertToImagej(image, title, true);
    }

    public static ImagePlus convertToImagej(ImageDouble image, String title, boolean useLogarithm) {
        if (image.isPSD()) {
            image.convertPSD(useLogarithm);
        }

        int w = image.getXsize();
        int h = image.getYsize();
        int d = image.getZsize();
        long n = image.getNsize();

        ImagePlus ip = convertToImagej(image.getData(), w, h, d, n, title);

        // Sets associated file info.
        File f = new File(title);
        FileInfo fi = new FileInfo();
        fi.directory = f.getParent();
        fi.fileName = f.getName();
        ip.setFileInfo(fi);

        return ip;
    }

    public static ImagePlus convertToImagej(Projection projection, String title) {

        int w = projection.getXsize();
        int h = projection.getYsize();
        int d = projection.getZsize();

        return convertToImagej(projection.getData(), w, h, d, 1, title);
    }

    private static ImagePlus convertToImagej(double array[], int w, int h, int d, long n, String title) {
        int sliceSize = w * h;
        int imageSize = sliceSize * d;
        double out[] = new double[sliceSize];
        ImageStack is = new ImageStack(w, h);

        for (int i = 0; i < n; i++) {
            int offset = i * imageSize;
            for (int j = 0; j < d; j++) {
                System.arraycopy(array, offset + j * sliceSize, out, 0, sliceSize);

                FloatProcessor processor = new FloatProcessor(w, h, out);
                is.addSlice(String.valueOf(i), processor);
            }
        }

        ImagePlus ip = new ImagePlus(title, is);

        // Normalize by default
        StackStatistics ss = new StackStatistics(ip);
        ip.getProcessor().setMinAndMax(ss.min, ss.max);

        return ip;
    }

    public static ImagePlus convertToImagej(MetaData md) {
        LinkedList<String> missing = new LinkedList<String>();
        ImagePlus ip = null;

        if (md.containsLabel(MDLabel.MDL_IMAGE)) {
            ImageStack is = null;

            long ids[] = md.findObjects();

            for (long id : ids) {
                String filename = md.getValueString(MDLabel.MDL_IMAGE, id, true);

                try {
                    ImageDouble img = new ImageDouble(filename);
                    ImagePlus slice = ImageConverter.convertToImagej(img, filename);

                    if (is == null) {
                        is = new ImageStack(slice.getWidth(), slice.getHeight());
                    }

                    is.addSlice(filename, slice.getProcessor());
                } catch (Exception ex) {
                    missing.add(filename);
                    //ex.printStackTrace();
                }
            }

            ip = new ImagePlus(md.getFilename(), is);

            // Sets associated file info.
            File f = new File(md.getFilename());
            FileInfo fi = new FileInfo();
            fi.directory = f.getParent();
            fi.fileName = f.getName();
            ip.setFileInfo(fi);

            // Normalize by default
            StackStatistics ss = new StackStatistics(ip);
            ip.getProcessor().setMinAndMax(ss.min, ss.max);
        }

        // Tells user about missing files.
        if (!missing.isEmpty()) {
            String message = "There are missing files:\n";
            for (int i = 0; i < missing.size(); i++) {
                message += missing.get(i) + "\n";
            }

            IJ.error(message);
        }

        return ip;
    }
    /*
    public static ImagePlus convertToImagePlus(Vector<TableImageItem> items, String path) {
    ImagePlus imp = convertToImagePlus(items);

    // Sets associated file info.
    File f = new File(path);
    FileInfo fi = new FileInfo();
    fi.directory = f.getParent();
    fi.fileName = f.getName();
    imp.setFileInfo(fi);

    return imp;
    }*/

    public static ImagePlus convertToImagePlus(Vector<AbstractTableImageItem> items) {
        ImageStack is = null;

        for (int i = 0; i < items.size(); i++) {
            AbstractTableImageItem item = items.elementAt(i);

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

    public static void revert(ImagePlus ip, String path) throws Exception {
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

        ip.setStack(is);
    }

    public static ImageDouble convertToXmipp(ImagePlus ip) {
        ImageDouble image = new ImageDouble();

        int w = ip.getWidth();
        int h = ip.getHeight();
        int d = ip.getStackSize();

        double data[] = new double[w * h * d];
        for (int i = 0; i < d; i++) {
            float slice[] = (float[]) ip.getStack().getProcessor(i + 1).getPixels();
            System.arraycopy(Tools.toDouble(slice), 0, data, i * w * h, w * h);
        }
        try {
            image.setData(w, h, d, data);
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        return image;
    }
    /*
    public static ImagePlus toImagePlus(Vector<TableImageItem> items, String title) {
    TableImageItem item = items.elementAt(0);

    ImageStack is = new ImageStack(item.getWidth(), item.getHeight());

    for (int i = 0; i < items.size(); i++) {
    is.addSlice(String.valueOf(i), items.elementAt(i).getPreview().getProcessor());
    }

    return new ImagePlus(title, is);
    }*/
}
