/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.listitems;

import browser.Cache;
import browser.LABELS;
import browser.imageitems.ImageConverter;
import browser.imageitems.ImageDimension;
import ij.IJ;
import ij.ImagePlus;
import java.io.File;
import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public class XmippImageItem extends AbstractImageItem {

    protected final static int FIRST_IMAGE = 0;
    protected final static int FIRST_SLICE = 0;
    public final static int MID_SLICE = -1;
    public int slice = MID_SLICE, nimage = FIRST_IMAGE;

    public XmippImageItem(File file, Cache cache) {
        super(file, cache);

        loadImageData();    // Loads size at start up.
    }

    protected void loadImageData() {
        try {
            ImageDouble image = new ImageDouble();
            image.readHeader(file.getAbsolutePath());

            dimension = new ImageDimension(image);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    public String getKey() {
        return super.getKey() + "-" + slice + "-" + nimage;
    }

    protected ImagePlus loadPreview(int w, int h) {
        return loadPreview(w, h, slice, nimage);
    }

    protected ImagePlus loadPreview(int w, int h, int slice, int nimage) {
        ImagePlus ip = null;

        try {
            String path = file.getAbsolutePath();
            //String fileName = file.getName();
            ImageDouble image = new ImageDouble();
            //System.out.println(" *** Loading preview: " + path + " / w=" + w + " / h=" + h + " / d=" + slice + " n=" + nimage);

            double factor = getFactor(dimension.width, dimension.height, w, h);

            int w_ = (int) Math.ceil(dimension.width / factor);
            int h_ = (int) Math.ceil(dimension.height / factor);

            image.readPreview(path, w_, h_, slice, nimage);
            ip = ImageConverter.convertToImagej(image, path);
        } catch (Exception e) {
            e.printStackTrace();
            IJ.error(e.getMessage());
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
            ImageDouble image = new ImageDouble();

            image.read(path, nimage);
            ip = ImageConverter.convertToImagej(image, path);
        } catch (Exception e) {
            e.printStackTrace();
            IJ.error(e.getMessage());
        }

        return ip;
    }

    public String getImageInfo() {
        loadImageData();

        return "<html>"
                + LABELS.LABEL_WIDTH + dimension.width + "<br>"
                + LABELS.LABEL_HEIGHT + dimension.height + "<br>"
                + (dimension.depth > 1 ? LABELS.LABEL_DEPTH + dimension.depth + "<br>" : "")
                + (dimension.nimages > 1 ? "<hr>" + LABELS.LABEL_NIMAGES + dimension.nimages : "")
                + "</html>";
    }
}
