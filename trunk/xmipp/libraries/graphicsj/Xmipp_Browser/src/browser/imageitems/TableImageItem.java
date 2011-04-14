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
    protected double scale = 1.0;

    public TableImageItem(File file, Cache cache) {
        this(file, cache, ImageDouble.FIRST_SLICE);
    }

    public TableImageItem(File file, Cache cache, int slice) {
        this(file, cache, slice, ImageDouble.FIRST_IMAGE);
    }

    public TableImageItem(File file, Cache cache, int slice, int nimage) {
        super(file, cache);

        this.nslice = slice;
        this.nimage = nimage;
    }

    @Override
    public String getKey() {
        return super.getKey() + "[" + scale + "]";
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

    public boolean isEnabled() {
        return enabled;
    }

    public boolean isSelected() {
        return selected;
    }

    public void setZoomScale(double scale) {
        this.scale = scale;
    }

    public int getThumbnailWidth() {
        return (int) ((double) super.getWidth() * scale);
    }

    public int getThumbnailHeight() {
        return (int) ((double) super.getHeight() * scale);
    }

    public ImagePlus getPreview() {
        return getPreview(getThumbnailWidth(), getThumbnailHeight());
    }

    @Override
    public String getLabel() {
        String sliceStr = isVolume() ? String.valueOf(nslice) : "";
        String nimageStr = isStack() ? String.valueOf(nimage) : "";
        String delim = isVolume() && isStack() ? "/" : "";

        String extra = isVolume() || isStack() ? "[" + sliceStr + delim + nimageStr + "]" : "";

        return super.getLabel() + extra;
    }

    @Override
    public String toString() {
        return file.getName() + "@" + nimage + "#" + nslice;
    }
}
