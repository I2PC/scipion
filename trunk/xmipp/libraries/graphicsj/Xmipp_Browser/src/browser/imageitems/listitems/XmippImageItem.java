/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.listitems;

import browser.Cache;
import browser.imageitems.ImageConverter;
import ij.IJ;
import ij.ImagePlus;
import java.io.File;
import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public class XmippImageItem extends AbstractImageItem {

    public int slice = ImageDouble.MID_SLICE, nimage = ImageDouble.FIRST_IMAGE;

    public XmippImageItem(File file, Cache cache) {
        super(file, cache);

        loadImageData();    // Loads size at start up.
    }

    @Override
    public String getKey() {
        return super.getKey() + "_" + slice + "_" + nimage;
    }

    protected ImagePlus loadPreview(int w, int h) {
        return loadPreview(w, h, slice, nimage);
    }

    protected ImagePlus loadPreview(int w, int h, int slice, int nimage) {
        System.out.println(" >>> Loading preview: " + getKey());
        ImagePlus ip = null;

        try {
            String path = file.getAbsolutePath();
            //String fileName = file.getName();
            ImageDouble image = new ImageDouble();
            //System.out.println(" *** Loading preview: " + path + " / w=" + w + " / h=" + h + " / d=" + slice + " n=" + nimage);

            double factor = getFactor(getWidth(), getHeight(), w, h);

            int w_ = (int) Math.ceil(getWidth() / factor);
            int h_ = (int) Math.ceil(getHeight() / factor);

            image.readPreview(path, w_, h_, slice, nimage);
            ip = ImageConverter.convertToImagej(image, path);
        } catch (Exception e) {
            System.err.println(" >>> Error loading preview: " + getKey());
            e.printStackTrace();
//            IJ.error(e.getMessage());
        }

        return ip;
    }

    public String getFileName() {
        return file.getName();
    }

    private static double getFactor(int W, int H, int w, int h) {
        double factor;

        if (W > H) {
            factor = (double) W / (double) w;
        } else {
            factor = (double) H / (double) h;
        }

        return factor;
    }

    public ImagePlus getImagePlus() {
        return getImagePlus(nimage);
    }

    public ImagePlus getImagePlus(int nimage) {
        ImagePlus ip = null;

        try {
            String path = file.getAbsolutePath();
            System.out.println(" *** Reading ImagePlus [from disk]: " + getFileName() + " #image: " + nimage);
            ImageDouble image = new ImageDouble();

            image.read(path, nimage);
            ip = ImageConverter.convertToImagej(image, path);
        } catch (Exception e) {
            e.printStackTrace();
            IJ.error(e.getMessage());
        }

        return ip;
    }
}
