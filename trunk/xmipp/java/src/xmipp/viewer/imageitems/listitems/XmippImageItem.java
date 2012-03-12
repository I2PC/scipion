/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.imageitems.listitems;

import xmipp.utils.Cache;
import xmipp.utils.DEBUG;
import ij.IJ;
import ij.ImagePlus;
import java.io.File;
import xmipp.jni.ImageGeneric;
import xmipp.ij.commons.XmippImageConverter;

/**
 *
 * @author Juanjo Vega
 */
public class XmippImageItem extends AbstractImageItem {

    public int nslice = ImageGeneric.MID_SLICE;
    public long nimage = ImageGeneric.FIRST_IMAGE;

    public XmippImageItem(File file, Cache cache) {
        super(file, cache);
    }

    @Override
    public String getKey() {
        return super.getKey() + "_" + nslice + "_" + nimage;
    }

    protected ImagePlus loadPreview(int w, int h) {
        return loadPreview(w, h, nslice, nimage);
    }

    protected ImagePlus loadPreview(int w, int h, int slice, long nimage) {
        ImagePlus ip = null;

        try {
            String path = file.getAbsolutePath();
            ImageGeneric image = new ImageGeneric(path);

            double factor = getFactor(getWidth(), getHeight(), w, h);

            int w_ = (int) Math.ceil(getWidth() / factor);
            int h_ = (int) Math.ceil(getHeight() / factor);

            image.read(w_, h_, slice, nimage);
            ip = XmippImageConverter.readToImagePlus(image);
        } catch (Exception ex) {
            System.out.println("ERROR: loading preview: " + getKey());
        }

        return ip;
    }

    public static double getFactor(int W, int H, int w, int h) {
        double factor;

        if (W > H) {
            factor = (double) W / (double) w;
        } else {
            factor = (double) H / (double) h;
        }

        return factor;
    }

    public ImagePlus getImagePlus() {
        return getImagePlus(nimage, nslice);
    }

    public ImagePlus getImagePlus(long nimage) {
        return getImagePlus(nimage, ImageGeneric.ALL_SLICES);
    }

    public ImagePlus getImagePlus(int nslice) {
        return getImagePlus(ImageGeneric.FIRST_IMAGE, nslice);
    }

    public ImagePlus getImagePlus(long nimage, int nslice) {
        ImagePlus imp = null;

        if (exists()) {
            try {
                DEBUG.printMessage(" *** Reading ImagePlus [from disk]: " + getKey());
                String path = getAbsoluteFileName();
                ImageGeneric image = new ImageGeneric(path);

                //image.readPreview(path, getWidth(), getHeight(), nslice, nimage);
                image.read(nimage);
                imp = XmippImageConverter.readToImagePlus(image);
                imp.setTitle(getLabel());
                //getPreview(getWidth(), getHeight());
            } catch (Exception ex) {
                IJ.error(ex.getMessage());
//                ex.printStackTrace();
            }
        }

        return imp;
    }
}
