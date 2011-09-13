/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.tableitems;

import browser.Cache;
import browser.DEBUG;
import browser.ICONS_MANAGER;
import browser.imageitems.ImageConverter;
import browser.imageitems.ImageDimension;
import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageStatistics;
import java.io.File;
import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public abstract class AbstractTableImageItem {

    protected Cache cache;
    protected ImageDimension dimension;
    protected ImageStatistics statistics;
    protected boolean selected;
    protected double scale = 1.0;

    public AbstractTableImageItem(Cache cache) {
        this.cache = cache;
    }

    // CompareTo for mergesort algorithm.
    public int compareToByLabel(AbstractTableImageItem item, int label) {
        String a = String.valueOf(getLabelValue(label));
        String b = String.valueOf(item.getLabelValue(label));

        return a.compareTo(b);
    }

    public abstract String getPath();

    public abstract long getNImage();

    public abstract int getNSlice();

    protected void loadImageData() {
        try {
            ImageDouble image = new ImageDouble();
            image.readHeader(getPath());

            dimension = new ImageDimension(image);
        } catch (Exception ex) {
            //throw new RuntimeException(ex);
            DEBUG.printMessage(ex.getMessage() + ": " + getPath() + " (image=" + getNImage() + ", slice=" + getNSlice() + ")");
        }
    }

    public boolean exists() {
        File file = new File(getPath());

        return file.exists();
    }

    public abstract void setEnabled(boolean enabled);

    public abstract boolean isEnabled();

    public void setSelected(boolean selected) {
        this.selected = selected;
    }
//
//    public String getAbsoluteFileName() {
//        return Filename.getFilename(get);
//    }

    public boolean isSelected() {
        return selected;
    }

    public boolean isStack() {
        return dimension != null && dimension.getNimages() > 1;
    }

    public boolean isVolume() {
        return dimension != null && dimension.getDepth() > 1;
    }

//    public double getZoomScale_() {
//        return scale;
//    }
    public void setZoomScale(double scale) {
        this.scale = scale;
    }

    public int getWidth() {
        return dimension != null ? dimension.getWidth() : ICONS_MANAGER.DEFAULT_PREVIEW_WIDTH;
    }

    public int getHeight() {
        return dimension != null ? dimension.getHeight() : ICONS_MANAGER.DEFAULT_PREVIEW_HEIGHT;
    }

    public String getKey() {
        return getPath() + ":i" + getNImage() + ":s" + getNSlice() + ":w" + getThumbnailWidth() + ":h" + getThumbnailHeight();
    }

    public int getThumbnailWidth() {
        return (int) (getWidth() * scale);
    }

    public int getThumbnailHeight() {
        return (int) (getHeight() * scale);
    }

    public abstract String getTooltipText();

    public abstract String getTitle();

    public abstract Object getLabelValue(int label);

    public ImagePlus getImagePlus() {
        ImagePlus ip = null;

        if (exists()) {
            try {
                System.out.println(" *** Reading ImagePlus [from disk]: " + getKey());
                ImageDouble image = new ImageDouble();

                image.readPreview(getPath(), getWidth(), getHeight(), getNSlice(), getNImage());
                ip = ImageConverter.convertToImagej(image, getTitle());

                ip.setTitle(getTitle());
            } catch (Exception ex) {
                IJ.error(ex.getMessage());
            }
        }

        return ip;
    }

    public ImagePlus getPreview() {
        return getPreview(getThumbnailWidth(), getThumbnailHeight());
    }

    public ImagePlus getPreview(int w, int h) {
        ImagePlus preview;

        if (getWidth() > 0 && getHeight() > 0) {
            // Tries to load from cache.
            preview = (ImagePlus) cache.get(getKey());

            // If not in cache.
            if (preview == null) {
//                System.out.println("Reading from disk: " + getKey());
                preview = loadPreview(w, h);

                if (preview != null) {
                    cache.put(getKey(), preview);
                }
            }

            // Preview might be loaded from cache if it has been already
            // referenced by another item, but statistics might be still null.
            if (statistics == null) {
                /*statistics = preview.getStatistics(
                ImageStatistics.MIN_MAX + ImageStatistics.MEAN + ImageStatistics.STD_DEV);*/
                //System.out.println("Loading statistics for: " + getLabelAsString());
                int moptions = ImageStatistics.MIN_MAX + ImageStatistics.MEAN + ImageStatistics.STD_DEV;
                statistics = ImageStatistics.getStatistics(preview.getProcessor(), moptions, preview.getCalibration());
            }
        } else {    // Null preview.
            preview = ICONS_MANAGER.MISSING_ITEM;
        }

        return preview;
    }

    protected ImagePlus loadPreview(int w, int h) {
//        System.out.println(" >>> Loading preview: " + getKey());
        ImagePlus ip = null;

        try {
            //String fileName = file.getName();
            ImageDouble image = new ImageDouble();
            //System.out.println(" *** Loading preview: " + path + " / w=" + w + " / h=" + h + " / d=" + nslice + " n=" + nimage);

            double factor = getFactor(getWidth(), getHeight(), w, h);

            int w_ = (int) Math.ceil(getWidth() / factor);
            int h_ = (int) Math.ceil(getHeight() / factor);

            image.readPreview(getPath(), w_, h_, getNSlice(), getNImage());
            ip = ImageConverter.convertToImagej(image, getTitle());
        } catch (Exception ex) {
            ip = ICONS_MANAGER.MISSING_ITEM;
        }

        return ip;
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

    public ImageStatistics getStatistics() {
        return statistics;
    }
}
