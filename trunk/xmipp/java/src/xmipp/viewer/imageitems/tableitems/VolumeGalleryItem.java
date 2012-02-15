/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.imageitems.tableitems;

import xmipp.utils.Cache;
import xmipp.utils.XmippLabel;
import java.io.File;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;

/**
 *
 * @author Juanjo Vega
 */
public class VolumeGalleryItem extends AbstractGalleryImageItem {

    protected String filename;
    protected int slice;
    protected boolean enabled = true;

    public VolumeGalleryItem(String filename, int slice, Cache cache) {
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
        return ImageGeneric.FIRST_IMAGE;
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
        return filename + " / slice: " + getNSlice();
    }

    @Override
    public String getTooltipText() {
        return getAbsoluteFileName() + " (" + slice + ")";
    }

    public Object getLabelValue(int label) {

        switch (label) {
            case MDLabel.MDL_IMAGE:
                int digits = String.valueOf(dimension.getDepth()).length();
                String sliceStr = String.format("%1$0" + digits + "d", getNSlice());

                return sliceStr + "#" + filename;
            case MDLabel.MDL_ENABLED:
                return enabled;
        }

        return XmippLabel.LABEL_UNKNOWN;
    }
}
