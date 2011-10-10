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
//        System.out.println(" >>> Loading preview: " + getKey());
        ImagePlus ip = null;

        try {
            String path = file.getAbsolutePath();
            //String fileName = file.getName();
            ImageDouble image = new ImageDouble();
            //System.out.println(" *** Loading preview: " + path + " / w=" + w + " / h=" + h + " / d=" + nslice + " n=" + nimage);

            double factor = XmippImageItem.getFactor(getWidth(), getHeight(), w, h);

            int w_ = (int) Math.ceil(getWidth() / factor);
            int h_ = (int) Math.ceil(getHeight() / factor);

            image.readPreview(path, w_, h_);
            ip = ImageConverter.convertToImageJ(image, path);
        } catch (Exception ex) {
            System.err.println(" >>> Error loading preview: " + getKey());
//            ex.printStackTrace();
            //throw new RuntimeException(ex);
        }

        return ip;
    }
}
