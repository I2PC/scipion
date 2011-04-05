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
    //protected double zoom = 1.0;
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

    public void setNormalized(double min, double max) {
        this.normalized = true;

        this.min = min;
        this.max = max;
    }

    public void resetNormalized() {
        this.normalized = false;
    }

    public boolean isEnabled() {
        return enabled;
    }

    public boolean isSelected() {
        return selected;
    }

    /*    public ImagePlus getPreview() {
    return getPreview(getThumbnailWidth(), getThumbnailHeight());
    }*/
    public ImagePlus getPreview() {
        return getPreview(getWidth(), getHeight());
    }

    public ImagePlus getPreview(double zoom) {
        ImagePlus preview = getPreview(
                (int) ((double) getWidth() * zoom),
                (int) ((double) getHeight() * zoom));

        // If normalization is active...
        if (normalized) {
            preview.getProcessor().setMinAndMax(min, max);
        }

        return preview;
    }

    @Override
    public String getLabel() {
        String sliceStr = isVolume() ? String.valueOf(slice) : "";
        String nimageStr = isStack() ? String.valueOf(nimage) : "";
        String delim = isVolume() && isStack() ? "/" : "";

        String extra = isVolume() || isStack() ? "[" + sliceStr + delim + nimageStr + "]" : "";

        return super.getLabel() + extra;
    }

    /*    public int getThumbnailWidth() {
    return (int) ((double) getWidth() * zoom);
    }

    public int getThumbnailHeight() {
    return (int) ((double) getHeight() * zoom);
    }*/
    public String getTooltipText() {
        /*        String sliceStr = isVolume() ? " slice: " + slice : "";
        String nimageStr = isStack() ? " image: " + nimage : "";

        return getLabel() + sliceStr + nimageStr;*/
        return getLabel();
    }

    @Override
    public String toString() {
        return file.getName() + "@" + nimage + "#" + slice;
    }
}
