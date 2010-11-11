/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems;

import browser.Cache;
import browser.imageitems.listitems.AbstractImageItem;
import browser.table.ImagesTableModel;
import ij.IJ;
import ij.ImagePlus;
import xmipp.io.ij.FileOpener;

/**
 *
 * @author Juanjo Vega
 */
public class TableImageItem {

    protected boolean enabled = true;
    protected boolean selected = false;
    protected double zoom = 1.0;
    protected boolean normalized = false;
    protected Cache cache;
    protected AbstractImageItem item;
    protected ImagesTableModel imagesTableModel;
    //protected double min, max;

    public TableImageItem(AbstractImageItem item, Cache cache, ImagesTableModel imagesTableModel) {
        this(item, cache, FileOpener.MID_SLICE, imagesTableModel);
    }

    public TableImageItem(AbstractImageItem item, Cache cache, int slice, ImagesTableModel imagesTableModel) {
//        super(item.getFile(), cache, slice);

        this.item = item;
        this.cache = cache;
        this.imagesTableModel = imagesTableModel;
    }

    public String getKey() {
        return item.getKey();
    }

    public String getLabel() {
        return item.getLabel();
    }

    public String getFileName() {
        return item.getFileName();
    }

    public void setEnabled(boolean enabled) {
        this.enabled = enabled;
        removeFromCache();  // Forces repainting
    }

    public void setSelected(boolean selected) {
        this.selected = selected;
    }

    public void setNormalized(boolean normalized/*double min, double max*/) {
        this.normalized = normalized;
//        this.min = min;
//        this.max = max;
    }

    public boolean isEnabled() {
        return enabled;
    }

    public boolean isSelected() {
        return selected;
    }

    protected void removeFromCache() {
        cache.remove(getKey());
    }

    public void setZoom(int zoom) {
        this.zoom = (double) zoom / 100.0;
    }

    public ImagePlus getPreview() {
        return getPreview(getThumbnailWidth(), getThumbnailHeight());
    }

    public ImagePlus getPreview(int w, int h) {
        ImagePlus preview = (ImagePlus) cache.get(getKey());

        // If not in cache...
        if (preview == null) {
            preview = item.loadPreview(w, h, item.slice);   // ...loads it from disk.

            // If item has been normalized, applies it.
            if (normalized) {
                preview.getProcessor().setMinAndMax(imagesTableModel.getMin(), imagesTableModel.getMax());
                preview.updateImage();
            }

            cache.put(getKey(), preview);   // Stores it
        }

        return preview;
    }

    public int getWidth() {
        return item.width;
    }

    public int getHeight() {
        return item.height;
    }

    public int getThumbnailWidth() {
        return (int) ((double) item.width * zoom);
    }

    public int getThumbnailHeight() {
        return (int) ((double) item.height * zoom);
    }

    public ImagePlus getImagePlus() {
        ImagePlus ip_ = IJ.openImage(item.getFile().getAbsolutePath());

        ImagePlus ip;
        if (item.slice != FileOpener.MID_SLICE) {
            ip = new ImagePlus();
            ip.setProcessor("", ip_.getStack().getProcessor(item.slice));
        } else {
            ip = ip_;
        }

        return ip;
    }

    public boolean loadImageData() {
        return item.loadImageData();
    }

    protected ImagePlus loadPreview(int w, int h, int slice) {
        return item.getPreview(w, h);
    }
}
