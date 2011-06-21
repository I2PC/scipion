/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.tableitems;

import browser.Cache;
import java.io.File;
import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public class VolumeTableItem extends AbstractTableImageItem {

    protected String filename;
    protected int slice;
    protected boolean enabled = true;

    public VolumeTableItem(String filename, int slice, Cache cache) {
        super(cache);

        this.filename = filename;
        this.slice = slice;

        loadImageData();
    }

    @Override
    public boolean isEnabled() {
        return enabled;
    }

    @Override
    public void setEnabled(boolean enabled) {
        this.enabled = enabled;
    }

    public String getAbsoluteFileName() {
        File f = new File(filename);

        return f.getAbsolutePath();
    }

    @Override
    public long getNImage() {
        return ImageDouble.FIRST_IMAGE;
    }

    @Override
    public int getNSlice() {
        return slice;
    }

    @Override
    public String getPath() {
        return filename;
    }

    @Override
    public String getTitle() {
        return filename + ": " + dimension.getDepth() + " slices.";
    }

    @Override
    public String getTooltipText() {
        return getAbsoluteFileName() + "[" + slice + "]";
    }

    public String getLabel() {
        return filename + "[" + slice + "]";
    }
}
