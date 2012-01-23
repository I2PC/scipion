/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.gallery.models;

import xmipp.utils.DEBUG;
import xmipp.viewer.imageitems.tableitems.AbstractGalleryImageItem;
import xmipp.viewer.imageitems.tableitems.VolumeGalleryItem;
import java.util.ArrayList;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class VolumeTableModel extends AbstractXmippTableModel {

    public VolumeTableModel() {
        super();
    }

    public VolumeTableModel(String filename) {
        super(filename);
    }

    protected String populateTable(String filename) {
        String message = null;

        try {
        	setCacheSize(filename);
            ImageGeneric image = new ImageGeneric(filename);

            int nslices = image.getZDim();
            long nimages = image.getNDim();

            for (long i = ImageGeneric.FIRST_IMAGE; i <= nimages; i++) {
                for (int j = ImageGeneric.FIRST_SLICE; j <= nslices; j++) {
                    data.add(new VolumeGalleryItem(filename, j, cache));
                }
            }
        } catch (Exception ex) {
            message = ex.getMessage();
        }

        return message;
    }

    @Override
    public String[] getLabels() {
        labelsValues = new int[]{MDLabel.MDL_IMAGE, MDLabel.MDL_ENABLED};
        String labels[] = new String[labelsValues.length];

        try {
            for (int i = 0; i < labelsValues.length; i++) {
                labels[i] = MetaData.label2Str(labelsValues[i]);
            }
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }

        return labels;
    }

    @Override
    protected void getMinAndMax() {
        try {
            ImageGeneric ig = new ImageGeneric(filename);
            ig.read(ImageGeneric.FIRST_IMAGE);
            double stats[] = ig.getStatistics();

            min = stats[0];
            max = stats[1];
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        DEBUG.printMessage(" *** VolumeTableModel: m=" + min + " M=" + max);
    }

    public String getFilename() {
        return filename;
    }

    public String getTitle() {
        String strImageSize = "";
        if (getSize() > 0) {
            AbstractGalleryImageItem item = getAllItems().get(0);
            strImageSize = " (" + item.getWidth() + " x " + item.getHeight() + ")";
        }

        return filename + ": " + getSize() + " slices." + strImageSize;
    }

    @Override
    public boolean isStack() {
        return false;
    }

    @Override
    public boolean isVolume() {
        return true;
    }

    @Override
    public boolean isMetaData() {
        return false;
    }

    @Override
    public boolean containsGeometryInfo() {
        return false;
    }

    @Override
    public boolean saveAsMetadata(String path, boolean all) {
        try {
            MetaData md = new MetaData();

            md.addLabel(MDLabel.MDL_ENABLED);
            md.addLabel(MDLabel.MDL_IMAGE);
            ArrayList<AbstractGalleryImageItem> items = all ? data : getSelectedItems();

            for (int i = 0; i < items.size(); i++) {
                AbstractGalleryImageItem item = items.get(i);

                long id = md.addObject();

                String image = (String) item.getLabelValue(MDLabel.MDL_IMAGE);
                int enabled = item.isEnabled() ? 1 : 0;

                md.setValueInt(MDLabel.MDL_ENABLED, enabled, id);
                md.setValueString(MDLabel.MDL_IMAGE, image, id);
            }

            md.write(path);

            return true;
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }

        return false;
    }
}
