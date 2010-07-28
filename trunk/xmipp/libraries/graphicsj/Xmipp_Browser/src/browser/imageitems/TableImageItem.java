/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems;

import browser.Cache;
import ij.IJ;
import ij.ImagePlus;
import java.io.File;
import xmipp.Spider_Reader;

/**
 *
 * @author Juanjo Vega
 */
public class TableImageItem extends ImageItem {

    protected boolean selected = false;
    protected double zoom = 1.0;
    protected boolean normalized = false;
    protected double min, max;

    public TableImageItem(File file, Cache cache) {
        super(file, cache);
    }

    public TableImageItem(File file, Cache cache, int slice) {
        super(file, cache, slice);
    }

    public String getKey() {
        return file.getAbsolutePath() + (slice >= 0 ? "-" + slice : "");
    }

    @Override
    public String getLabel() {
        String label = super.getLabel();

        if (slice >= 0) {
            label += " (slice " + slice + ")";
        }

        return label;
    }

    public void setNormalized(double min, double max) {
        normalized = true;
        this.min = min;
        this.max = max;
    }

    public void resetNormalized() {
        normalized = false;
    }

    public void setZoom(int zoom) {
        this.zoom = (double) zoom / 100.0;
    }

    public ImagePlus getPreview() {
        return getPreview(getThumbnailWidth(), getThumbnailHeight());
    }

    @Override
    public ImagePlus getPreview(int w, int h) {
        ImagePlus preview = null;

        if (w > 1 && h > 1) {
            boolean fromdisk = false;

            // If not in cache...
            if (cache.get(getKey()) == null) {
                fromdisk = true;
            }

            // ...it will be loaded from disk by super class, so...
            preview = super.getPreview(w, h);

            // ... if item has been normalized, then applies it when reloading from disk.
            if (fromdisk && normalized) {
                preview.getProcessor().setMinAndMax(min, max);
                preview.updateImage();
            }
        }

        return preview;
    }

    public void setSelected(boolean selected) {
        this.selected = selected;
    }

    public boolean isSelected() {
        return selected;
    }

    public int getThumbnailWidth() {
        return (int) ((double) width * zoom);
    }

    public int getThumbnailHeight() {
        return (int) ((double) height * zoom);
    }

    public ImagePlus getImagePlus() {
        ImagePlus ip_ = IJ.openImage(getFile().getAbsolutePath());//ImagesWindowFactory.openImageWindow(getFile());

        ImagePlus ip;
        if (slice != Spider_Reader.MID_SLICE) {
            ip = new ImagePlus();
            ip.setProcessor("", ip_.getStack().getProcessor(slice));
        } else {
            ip = ip_;
        }

        return ip;
    }
}
