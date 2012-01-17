/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.gallery.models;

import xmipp.utils.Cache;
import xmipp.utils.DEBUG;
import xmipp.utils.LABELS;
import xmipp.viewer.imageitems.tableitems.AbstractGalleryImageItem;
import xmipp.viewer.imageitems.tableitems.GalleryImageItem;
import ij.IJ;
import java.util.ArrayList;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class MDTableModel extends AbstractXmippTableModel {

    MetaData md;
    boolean containsGeometryInfo;

    public MDTableModel(String filename) {
        super(filename);
    }

    /*
     * For an array of single images.
     */
    public MDTableModel(String filenames[]) {
        super();

        String message = null;
        try {
            md = new MetaData();
            md.addLabel(MDLabel.MDL_IMAGE);
            md.addLabel(MDLabel.MDL_ENABLED);

            for (int i = 0; i < filenames.length; i++) {
                long id = md.addObject();

                md.setValueString(MDLabel.MDL_IMAGE, filenames[i], id);
                md.setValueInt(MDLabel.MDL_ENABLED, 1, id);
            }
        } catch (Exception ex) {
            message = ex.getMessage();
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

    void setCacheSize(MetaData md) throws Exception {
        // Calculates cache elements size.
        String firstImage = md.getValueString(MDLabel.MDL_IMAGE, md.firstObject(), true);
        ImageGeneric image = new ImageGeneric(firstImage);

        int imageSize = image.getXDim() * image.getYDim() * Cache.MAXPXSIZE;
        int elements = Cache.MEMORY_SIZE / imageSize;

        cache.resize(elements > 0 ? elements : 1);
    }

    @Override
    protected String populateTable(String path) {
        String message = null;

        try {
            if (md == null) {
                md = new MetaData(path);
            }

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
        long ids[] = null;

        try {
            setCacheSize(md);

            ids = md.findObjects();

            for (long id : ids) {
                GalleryImageItem item = new GalleryImageItem(id, md, MDLabel.MDL_IMAGE, cache);
                data.add(item);
            }
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }
    }

    public void setUseGeometry(boolean use) {
        for (int i = 0; i < getSize(); i++) {
            ((GalleryImageItem) data.get(i)).setReadGeo(use);
        }
        cache.clear();
    }

    @Override
    public String[] getLabels() {
        String labels[] = null;

        try {
            labelsValues = md.getActiveLabels();
            labels = new String[labelsValues.length];

            for (int i = 0; i < labelsValues.length; i++) {
                labels[i] = MetaData.label2Str(labelsValues[i]);
            }
        } catch (Exception ex) {
            DEBUG.printException(ex);
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
            AbstractGalleryImageItem item = getAllItems().get(0);
            strImageSize = " (" + item.getWidth() + " x " + item.getHeight() + ")";
        }

        String title = filename;//md.getFilename();
        return (title == null ? LABELS.TITLE_UNTITLED : title) + ": " + getSize() + " images." + strImageSize;
    }

    @Override
    protected void getMinAndMax() {
        try {
            double stats[] = md.getStatistics(false);

            min = stats[0];
            max = stats[1];
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }
    }

    @Override
    public boolean isStack() {
        boolean isStack = false;

        try {
            isStack = Filename.isStack(getFilename());
        } catch (Exception ex) {
        }

        return isStack;
    }

    @Override
    public boolean isVolume() {
        return false;
    }

    @Override
    public boolean isMetaData() {
        return Filename.isMetadata(getFilename());
    }

    @Override
    public boolean containsGeometryInfo() {
        return containsGeometryInfo;
    }

    @Override
    public boolean saveAsMetadata(String path, boolean all) {
        try {
            // Copies items into the new MetaData.
            ArrayList<AbstractGalleryImageItem> items = all ? data : getSelectedItems();
            long ids[] = new long[items.size()];

            for (int i = 0; i < items.size(); i++) {
                long id = ((GalleryImageItem) items.get(i)).getID();
                ids[i] = id;
            }

            MetaData output = new MetaData();
            output.importObjects(md, ids);

            output.write(path);

            return true;
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }

        return false;
    }
}
