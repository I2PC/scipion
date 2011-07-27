/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.tableitems;

import xmipp.Filename;
import xmipp.MDLabel;
import xmipp.MetaData;
import browser.Cache;
import browser.LABELS;
import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public class MDTableItem extends AbstractTableImageItem {

    protected MetaData md;
    protected long id;
    protected String originalValue;
    protected String path;
    protected long nimage;

    public MDTableItem(long id, MetaData md, Cache cache) {
        super(cache);

        this.id = id;
        this.md = md;

        originalValue = md.getValueString(MDLabel.MDL_IMAGE, id);
        String field = md.getValueString(MDLabel.MDL_IMAGE, id, true);
        this.path = Filename.getFilename(field);
        nimage = Filename.getNimage(field);

        loadImageData();
    }

    public boolean isEnabled() {
        return getValueInt(MDLabel.MDL_ENABLED) == 1;
    }

    public void setEnabled(boolean enabled) {
        md.setValueInt(MDLabel.MDL_ENABLED, enabled ? 1 : 0, id);
    }

    public String getAbsoluteFileName() {
        return Filename.getFilename(path);
    }

    public long getNImage() {
        return nimage;
    }

    @Override
    public String getPath() {
        return path;
    }

    @Override
    public int getNSlice() {
        return ImageDouble.FIRST_SLICE;
    }

    @Override
    public String getTitle() {
        return originalValue;
    }

    public String getTooltipText() {
        return path + (Filename.isStack(path) ? " [" + getNImage() + "]" : "");
    }

    public Object getLabelValue(int label) {
        Class class_ = MetaData.getLabelType(label);

        // Special label.
        if (label == MDLabel.MDL_ENABLED) {
            return md.getValueInt(label, id) > 0 ? Boolean.TRUE : Boolean.FALSE;
        }

        // Rest of them...
        if (class_ == Integer.class) {
            return md.getValueInt(label, id);
        }
        if (class_ == Double.class) {
            return md.getValueDouble(label, id);
        }
        if (class_ == String.class) {
            return md.getValueString(label, id);
        }

        return LABELS.LABEL_UNKNOWN;
    }

    // Metadata related methods.
    protected String getValueString(int label) {
        return md.getValueString(label, id, true);
    }

    protected boolean getValueBoolean(int label) {
        return md.getValueBoolean(label, id);
    }

    protected double getValueDouble(int label) {
        return md.getValueDouble(label, id);
    }

    protected double getValueInt(int label) {
        return md.getValueInt(label, id);
    }
}
