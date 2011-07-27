/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package metadata.images;

import browser.Cache;
import browser.imageitems.tableitems.AbstractTableImageItem;
import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public class TableImageItem extends AbstractTableImageItem {

    protected boolean enabled = true;
    protected String path;
    protected String originalValue;

    public TableImageItem(String path, String originalValue, Cache cache) {
        super(cache);

        this.path = path;
        this.originalValue = originalValue;

        loadImageData();
    }

    public void setEnabled(boolean enabled) {
        this.enabled = enabled;
    }

    public boolean isEnabled() {
        return enabled;
    }

    @Override
    public String getPath() {
        return path;
    }

    public void setPath(String path) {
        this.path = path;
    }

    public String getOriginalValue() {
        return originalValue;
    }

    public void setOriginalValue(String value) {
        originalValue = value;
    }
    @Override
    public long getNImage() {
        return ImageDouble.FIRST_IMAGE;
    }

    @Override
    public int getNSlice() {
        return ImageDouble.FIRST_SLICE;
    }

    @Override
    public String getTooltipText() {
        return path;
    }

    @Override
    public String getTitle() {
        return originalValue;
    }

    @Override
    public Object getLabelValue(int label) {
        return getOriginalValue();
    }

    @Override
    public String toString() {
        return getPath();
    }
}
