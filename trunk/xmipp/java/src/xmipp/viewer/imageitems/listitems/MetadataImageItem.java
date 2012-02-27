/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.imageitems.listitems;

import xmipp.utils.Cache;
import ij.IJ;
import ij.ImagePlus;
import java.io.File;
import xmipp.jni.ImageGeneric;
import xmipp.ij.XmippImageConverter;

/**
 *
 * @author Juanjo Vega
 */
public class MetadataImageItem extends AbstractImageItem {

    String originalFilename;

    public MetadataImageItem(File file, String originalFilename, Cache cache) {
        super(file, cache);

        this.originalFilename = originalFilename;
    }

    public String getValue() {
        return originalFilename;
    }

    @Override
    public String toString() {
        return originalFilename;
    }

    @Override
    public ImagePlus getImagePlus() {
        return IJ.openImage(file.getAbsolutePath());
    }

    @Override
    protected ImagePlus loadPreview(int w, int h) {
        ImagePlus ip = null;

        try {
            String path = file.getAbsolutePath();
            //String fileName = file.getName();
            ImageGeneric image = new ImageGeneric(path);

            double factor = XmippImageItem.getFactor(getWidth(), getHeight(), w, h);

            int w_ = (int) Math.ceil(getWidth() / factor);
            int h_ = (int) Math.ceil(getHeight() / factor);

            image.read( w_, h_);
            ip = XmippImageConverter.readImageGenericToImageJ(image);
        } catch (Exception ex) {
            System.err.println(" >>> Error loading preview: " + getKey());
//            ex.printStackTrace();
            //throw new RuntimeException(ex);
        }

        return ip;
    }
}
