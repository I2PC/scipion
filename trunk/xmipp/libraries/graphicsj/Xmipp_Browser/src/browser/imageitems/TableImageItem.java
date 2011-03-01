/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems;

import browser.Cache;
import browser.imageitems.listitems.XmippImageItem;
import ij.ImagePlus;
import java.io.File;
import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public class TableImageItem extends XmippImageItem {

    protected boolean enabled = true;
    protected boolean selected = false;
    protected double zoom = 1.0;
    protected boolean normalized = false;
    //protected ImagesTableModel imagesTableModel;
    //protected double min, max;

    public TableImageItem(File file, Cache cache) {//, ImagesTableModel imagesTableModel
        this(file, cache, ImageDouble.FIRST_SLICE);//, imagesTableModel);
    }

    public TableImageItem(File file, Cache cache, int slice) {//, ImagesTableModel imagesTableModel) {
        this(file, cache, slice, ImageDouble.FIRST_IMAGE);//, imagesTableModel);
    }

    public TableImageItem(File file, Cache cache, int slice, int nimage) {//, ImagesTableModel imagesTableModel) {
        super(file, cache);

        this.slice = slice;
        this.nimage = nimage;
//        this.imagesTableModel = imagesTableModel;
    }

    @Override
    public String getFileName() {
        return file.getAbsolutePath();
    }

    public void setEnabled(boolean enabled) {
        this.enabled = enabled;
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

    public void setZoom(int zoom) {
        this.zoom = (double) zoom / 100.0;
    }

    public ImagePlus getPreview() {
        return getPreview(getThumbnailWidth(), getThumbnailHeight());
    }

    public int getThumbnailWidth() {
        return (int) ((double) getWidth() * zoom);
    }

    public int getThumbnailHeight() {
        return (int) ((double) getHeight() * zoom);
    }

    public String getTooltipText() {
        String sliceStr = isVolume() ? " slice: " + slice : "";
        String nimageStr = isStack() ? " image: " + nimage : "";

        return getLabel() + sliceStr + nimageStr;
    }
}
