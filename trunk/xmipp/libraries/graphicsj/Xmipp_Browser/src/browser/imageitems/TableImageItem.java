/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems;

import browser.Cache;
import browser.table.ImagesTableModel;
import ij.IJ;
import ij.ImagePlus;
import java.io.File;
import xmipp.Spider_Reader;

/**
 *
 * @author Juanjo Vega
 */
public class TableImageItem extends ImageItem {

    public boolean enabled = true;
    public boolean selected = false;
    protected double zoom = 1.0;
    protected boolean normalized = false;
    protected ImagesTableModel imagesTableModel;
    //protected double min, max;

    public TableImageItem(File file, Cache cache, ImagesTableModel imagesTableModel) {
        //super(file, cache);
        this(file, cache, Spider_Reader.MID_SLICE, imagesTableModel);
    }

    public TableImageItem(File file, Cache cache, int slice, ImagesTableModel imagesTableModel) {
        super(file, cache, slice);
        this.imagesTableModel = imagesTableModel;
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

    public void setNormalized(boolean normalized/*double min, double max*/) {
        this.normalized = normalized;
//        this.min = min;
//        this.max = max;
    }
    /*
    public void resetNormalized() {
    normalized = false;
    }*/

    public void setZoom(int zoom) {
        this.zoom = (double) zoom / 100.0;
    }

    public ImagePlus getPreview() {
        return getPreview(getThumbnailWidth(), getThumbnailHeight());
    }

    @Override
    public ImagePlus getPreview(int w, int h) {
        ImagePlus preview = null;

        boolean fromdisk = false;

        // If not in cache...
        if (cache.get(getKey()) == null) {
            fromdisk = true;
        }

        // ...it will be loaded from disk by super class, so...
        preview = getZoomedPreview(w, h);

        // ... if item has been normalized, then applies it when reloading from disk.
        if (fromdisk && normalized) {
            preview.getProcessor().setMinAndMax(imagesTableModel.getMin(), imagesTableModel.getMax());
            preview.updateImage();
        }

        return preview;
    }

    private ImagePlus getZoomedPreview(int w, int h) {
        double tmp_zoom = zoom; // Stores zoom temporally.

        if (zoom > 1.0) {   // Max zoom to read from disk is 100.
            zoom = 1.0;
        }

        ImagePlus ip = super.getPreview(getThumbnailWidth(), getThumbnailHeight());

        // Zooms using ImageJ
        if (tmp_zoom > 1.0) {
            zoom = tmp_zoom;    // Restorers zoom.

            // Creates a new ImageProcessor by scaling the image, and replaces the original one.
            ip.setProcessor(ip.getTitle(), ip.getProcessor().resize(w, h));
        }

        return ip;
    }

    public int getThumbnailWidth() {
        return (int) ((double) width * zoom);
    }

    public int getThumbnailHeight() {
        return (int) ((double) height * zoom);
    }

    public ImagePlus getImagePlus() {
        ImagePlus ip_ = IJ.openImage(getFile().getAbsolutePath());

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
