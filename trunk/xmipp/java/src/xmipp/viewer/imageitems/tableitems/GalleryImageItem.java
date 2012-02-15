/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.imageitems.tableitems;

import xmipp.jni.Filename;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.Cache;
import xmipp.utils.DEBUG;
import xmipp.utils.XmippResource;
import xmipp.utils.XmippLabel;
import ij.ImagePlus;
import xmipp.jni.ImageGeneric;
import xmipp.ij.XmippImageConverter;

/**
 *
 * @author Juanjo Vega
 */
public class GalleryImageItem extends AbstractGalleryImageItem {

    protected long id;
    protected MetaData md;
    protected int label;
    protected String originalValue;
    protected String path;
    protected long nimage;
    protected boolean readGeo;

    public GalleryImageItem(long id, MetaData md, int label, Cache cache) {
        super(cache);

        this.id = id;
        this.md = md;
        this.label = label;

        try {
            originalValue = md.getValueString(label, id);
            if (originalValue != null) {
                String field = md.getValueString(label, id, true);
                path = Filename.getFilename(field);
                nimage = Filename.getNimage(field);

                loadImageData();
            }
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }
    }

    @Override
    protected ImagePlus loadPreview(int w, int h) {
        ImagePlus imp = XmippResource.MISSING_ITEM;

        try {
            if (readGeo && md.containsGeometryInfo()) {
                try {
                    ImageGeneric image = new ImageGeneric();
                    image.readApplyGeo(getAbsoluteFileName(), md, id, w, h);

                    imp = XmippImageConverter.convertImageGenericToImageJ(image);
                    imp.setTitle(getOriginalValue());
                } catch (Exception ex) {
                    DEBUG.printException(ex);
                    System.out.println("ERROR: readApplyGeo: " + getAbsoluteFileName() + " > " + ex.getMessage());
                }
            } else {
                imp = super.loadPreview(w, h);
            }
        } catch (Exception ex) {
            DEBUG.printException(ex);
            DEBUG.printMessage("fn: " + getAbsoluteFileName());
        }

        return imp;
    }

    public boolean isEnabled() {
        boolean enabled = true;

        try {
            enabled = getValueInt(MDLabel.MDL_ENABLED) == 1;
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }

        return enabled;
    }

    public void setEnabled(boolean enabled) {
        try {
            md.setValueInt(MDLabel.MDL_ENABLED, enabled ? 1 : 0, id);
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }
    }

    public void setReadGeo(boolean readGeo) {
        this.readGeo = readGeo;
    }

    public long getID() {
        return id;
    }

    public int getLabel() {
        return label;
    }

    public String getAbsoluteFileName() {
    	if (nimage != ImageGeneric.ALL_IMAGES)
    		return nimage + Filename.SEPARATOR + path;
    	return path;
    }

    public long getNImage() {
        return nimage;
    }

    @Override
    public String getPath() {
        return path;
    }

    public void setPath(String path) {
        this.path = path;
    }

    @Override
    public int getNSlice() {
        return ImageGeneric.FIRST_SLICE;
    }

    @Override
    public String getTitle() {
        return originalValue;
    }

    public String getOriginalValue() {
        return originalValue;
    }

    public void setOriginalValue(String value) {
        originalValue = value;
    }

    public String getTooltipText() {
        return originalValue;//path + (isStack() ? " [" + getNImage() + "]" : "");
    }

    public Object getLabelValue(int label) {
        try {
            Class class_ = MetaData.getLabelClass(label);

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
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }

        return XmippLabel.LABEL_UNKNOWN;
    }

    // Metadata related methods.
    protected String getValueString(int label) throws Exception {
        return md.getValueString(label, id, true);
    }

    protected boolean getValueBoolean(int label) throws Exception {
        return md.getValueBoolean(label, id);
    }

    protected double getValueDouble(int label) throws Exception {
        return md.getValueDouble(label, id);
    }

    protected double getValueInt(int label) throws Exception {
        return md.getValueInt(label, id);
    }

    @Override
    public String toString() {
        return originalValue;
    }
}
