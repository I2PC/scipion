/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.gallery.models;

import browser.DEBUG;
import browser.imageitems.tableitems.AbstractGalleryImageItem;
import browser.imageitems.tableitems.VolumeGalleryItem;
import java.util.ArrayList;
import xmipp.ImageDouble;
import xmipp.ImageGeneric;
import xmipp.MDLabel;
import xmipp.MetaData;

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
        String message = "";

        try {
            ImageDouble image = new ImageDouble();

            image.readHeader(filename);

            int nslices = image.getZsize();
            long nimages = image.getNsize();

            for (int i = ImageDouble.FIRST_IMAGE; i <= nimages; i++) {
                for (int j = ImageDouble.FIRST_SLICE; j <= nslices; j++) {
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

        for (int i = 0; i < labelsValues.length; i++) {
            labels[i] = MetaData.label2Str(labelsValues[i]);
        }

        return labels;
    }

    @Override
    protected void getMinAndMax() {
        try {
            ImageGeneric ig = new ImageGeneric(filename);
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
        MetaData md = new MetaData();

        md.addLabel(MDLabel.MDL_ENABLED);
        md.addLabel(MDLabel.MDL_IMAGE);

        try {
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
        }

        return false;
    }
}
