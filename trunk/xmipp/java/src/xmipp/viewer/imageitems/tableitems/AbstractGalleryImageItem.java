/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.imageitems.tableitems;

import xmipp.utils.Cache;
import xmipp.utils.DEBUG;
import xmipp.utils.Resources;
import xmipp.viewer.imageitems.ImageDimension;
import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageStatistics;
import java.io.File;
import xmipp.jni.ImageGeneric;
import xmipp.ij.XmippImageConverter;

/**
 *
 * @author Juanjo Vega
 */
public abstract class AbstractGalleryImageItem {

    protected Cache cache;
    protected ImageDimension dimension;
    protected ImageStatistics statistics;
    protected boolean selected;
    protected double scale = 1.0;

    public AbstractGalleryImageItem(Cache cache) {
        this.cache = cache;
    }

    // CompareTo for mergesort algorithm.
    public int compareToByLabel(AbstractGalleryImageItem item, int label) {
        String a = String.valueOf(getLabelValue(label));
        String b = String.valueOf(item.getLabelValue(label));

        return a.compareTo(b);
    }

    public abstract String getPath();

    public abstract long getNImage();

    public abstract int getNSlice();

    protected void loadImageData() {
        try {
            ImageGeneric image = new ImageGeneric(getPath());
            dimension = new ImageDimension(image);
            image.destroy();
        } catch (Exception ex) {
            //throw new RuntimeException(ex);
            DEBUG.printMessage(ex.getMessage() + ": " + getPath() + " (image=" + getNImage() + ", slice=" + getNSlice() + ")");
        }
    }

    public boolean exists() {
        String path = getPath();

        return path != null ? new File(path).exists() : false;
    }

    public abstract void setEnabled(boolean enabled);

    public abstract boolean isEnabled();

    public boolean isBiggerThan(int MAX_SIZE) {
        return getWidth() > MAX_SIZE || getHeight() > MAX_SIZE;
    }

    public void setSelected(boolean selected) {
        this.selected = selected;
    }

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
        return dimension != null ? dimension.getWidth() : Resources.DEFAULT_PREVIEW_WIDTH;
    }

    public int getHeight() {
        return dimension != null ? dimension.getHeight() : Resources.DEFAULT_PREVIEW_HEIGHT;
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
                //DEBUG.printMessage(" *** Reading ImagePlus [from disk]: " + getKey());
                ImageGeneric image = new ImageGeneric(getPath());

                ip = XmippImageConverter.convertToImageJ(image, getWidth(), getHeight(), getNSlice(), getNImage());
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
        String key = getKey();
        DEBUG.printMessage("getPreview: KEY: " + key);
        
        if (getWidth() > 0 && getHeight() > 0) {
            // Tries to load from cache.
            preview = (ImagePlus) cache.get(key);

            // If not in cache.
            if (preview == null) {
            	DEBUG.printMessage("   NOT IN CACHE, LOADING....KEY:" + key);
                preview = loadPreview(w, h);

                if (preview != null) {
                    cache.put(key, preview);
                    DEBUG.printMessage("   PUTTING IN CACHE....KEY:" + key);
                }
            }

            // Preview might be loaded from cache if it has been already
            // referenced by another item, but statistics might be still null.
            if (statistics == null) {
                int moptions = ImageStatistics.MIN_MAX + ImageStatistics.MEAN + ImageStatistics.STD_DEV;
                statistics = ImageStatistics.getStatistics(preview.getProcessor(), moptions, preview.getCalibration());
            }
        } else {    // Null preview.
            preview = Resources.MISSING_ITEM;
        }

        return preview;
    }

    protected ImagePlus loadPreview(int w, int h) {
        ImagePlus ip = Resources.MISSING_ITEM;
        String path = getPath();

        if (path != null) {
            try {
                double factor = getFactor(getWidth(), getHeight(), w, h);

                int w_ = (int) Math.ceil(getWidth() / factor);
                int h_ = (int) Math.ceil(getHeight() / factor);

                DEBUG.printMessage(String.format(" *** path: %s w=%d h=%d factor=%f s=%d, n=%d", 
                		getPath(), w_ , h_ , factor, getNSlice(), getNImage()));

                ImageGeneric image = new ImageGeneric(path);
                ip = XmippImageConverter.convertToImageJ(image, w_, h_, getNSlice(), getNImage());//,w_, h_, getNSlice(), getNImage());
                DEBUG.printMessage(String.format(" ***    size: %d", ip.getImageStackSize()));
                image.destroy();
            } catch (Exception ex) {
            	//DEBUG.printMessage("================== EXCEPTION ON GALLERYITEM loadPreview ===============");
                DEBUG.printException(ex);
            }
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
