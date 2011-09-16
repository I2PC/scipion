/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.models;

import browser.LABELS;
import browser.imageitems.tableitems.AbstractTableImageItem;
import browser.imageitems.tableitems.TableImageItem;
import ij.IJ;
import java.io.File;
import java.util.ArrayList;
import xmipp.MDLabel;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class MDTableModel extends AbstractXmippTableModel {

    MetaData md;
    boolean containsGeometryInfo = false;

    public MDTableModel(String filename) {
        super(filename);
    }

    /*
     * For an array of single images.
     */
    public MDTableModel(String filenames[]) {
        super();

        String message = null;

        for (int i = 0; i < filenames.length; i++) {
            String currentFile = filenames[i];

            File f = new File(currentFile);
            if (f.exists()) {
                try {
                    if (md == null) {
                        md = new MetaData();
                    }

                    long id = md.addObject();
                    md.setValueString(MDLabel.MDL_IMAGE, currentFile, id);
                } catch (Exception ex) {
                    message = ex.getMessage();
                }
            } else {
                message += "File not found: " + currentFile + "\n";
            }
        }

        populateTable();

        if (message != null) {
            IJ.error(message);
        }
    }

    public MDTableModel(String filenames[], boolean enabled[]) {
        this(filenames);

        // Set enabled/disabled
        if (enabled != null) {
            for (int i = 0; i < enabled.length; i++) {
                getAllItems().get(i).setEnabled(enabled[i]);
            }
        }
    }

    @Override
    protected String populateTable(String path) {
        String message = "";

        try {
            md = new MetaData(path);

            containsGeometryInfo = md.containsGeometryInfo();

            // Adds enabled field if not present.
            if (!md.containsLabel(MDLabel.MDL_ENABLED)) {
                md.addLabel(MDLabel.MDL_ENABLED);

                long ids[] = md.findObjects();
                for (long id : ids) {
                    md.setValueInt(MDLabel.MDL_ENABLED, 1, id);
                }
            }

            populateTable();
        } catch (Exception ex) {
            message = ex.getMessage();
        }

        return message;
    }

    private void populateTable() {
        // Populate table.
        long ids[] = md.findObjects();

        for (long id : ids) {
            data.add(new TableImageItem(id, md, MDLabel.MDL_IMAGE, cache));
        }
    }

    public void setUseGeometry(boolean use) {
        for (int i = 0; i < getSize(); i++) {
            ((TableImageItem) data.get(i)).setReadGeo(use);
        }
        cache.clear();
    }

    @Override
    public String[] getLabels() {
        labelsValues = md.getActiveLabels();
        String labels[] = new String[labelsValues.length];

        for (int i = 0; i < labelsValues.length; i++) {
            labels[i] = MetaData.label2Str(labelsValues[i]);
        }

        return labels;
    }

    @Override
    public String getFilename() {
        return md.getPath();
    }

    @Override
    public String getTitle() {
        String strImageSize = "";
        if (getSize() > 0) {
            AbstractTableImageItem item = getAllItems().get(0);
            strImageSize = " (" + item.getWidth() + " x " + item.getHeight() + ")";
        }

        String title = filename;//md.getFilename();
        return (title == null ? LABELS.TITLE_UNTITLED : title) + ": " + getSize() + " images." + strImageSize;
    }

    @Override
    protected void getMinAndMax() {
        double stats[] = md.getStatistics(false);

        min = stats[0];
        max = stats[1];
    }

    @Override
    public boolean isStack() {
        return true;
    }

    @Override
    public boolean isVolume() {
        return false;
    }

    @Override
    public boolean isMetaData() {
        return true;
    }

    @Override
    public boolean containsGeometryInfo() {
        return containsGeometryInfo;
    }

    @Override
    public boolean saveAsMetadata(String path, boolean all) {
        try {
            ArrayList<AbstractTableImageItem> items = getSelectedItems();
            for (int i = 0; i < items.size(); i++) {
                TableImageItem item = (TableImageItem) items.get(i);
                item.getID();
            }
            
            // md.copyRow(md2, id) ?
            
            md.write(path);

            return true;
        } catch (Exception ex) {
        }

        return false;
    }
}
