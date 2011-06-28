/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.models;

import browser.DEBUG;
import browser.imageitems.tableitems.VolumeTableItem;
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
        try {
            ImageDouble image = new ImageDouble();

            image.readHeader(filename);

            int nslices = image.getZsize();
            long nimages = image.getNsize();

            for (int i = ImageDouble.FIRST_IMAGE; i <= nimages; i++) {
                for (int j = ImageDouble.FIRST_SLICE; j <= nslices; j++) {
                    data.add(new VolumeTableItem(filename, j, cache));
                }
            }
        } catch (Exception ex) {
            return ex.getMessage();
        }

        return null;
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
        return filename + ": " + getSize() + " slices.";
    }

    @Override
    public boolean isStack() {
        return false;
    }

    @Override
    public boolean isVolume() {
        return true;
    }
}
