/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table;

import browser.DEBUG;
import browser.Cache;
import browser.imageitems.TableImageItem;
import ij.IJ;
import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;
import javax.swing.table.AbstractTableModel;
import xmipp.Filename;
import xmipp.ImageDouble;
import xmipp.MDLabel;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesTableModel extends AbstractTableModel {

    private String filename;
    private MetaData md;
    private long nimages = ImageDouble.FIRST_IMAGE;
    private int nslices = ImageDouble.FIRST_SLICE;
    //private double zoomScale = 1.0; // For zoom.
    private Vector<TableImageItem> data = new Vector<TableImageItem>();
    //private Vector<TableImageItem> selectedItems = new Vector<TableImageItem>();
    private LinkedList<TableImageItem> selectedItems = new LinkedList<TableImageItem>();
    private int rows, cols;
    private Cache cache = new Cache();
    private boolean normalize = false;
    private double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;

    public ImagesTableModel(String filename) {
        super();

        this.filename = filename;
        String message = null;

        File f = new File(filename);
        if (f.exists()) {
            if (Filename.isMetadata(filename)) {
                message = loadMetaData(filename);
            } else if (Filename.isSingleImage(filename)) {
                loadImages(new String[]{filename});
            } else if (Filename.isStackOrVolume(filename)) {
                message = loadStackOrVolume(filename);
            }

            if (message != null) {
                IJ.error(filename);
            }
        } else {
            IJ.error("File not found: " + filename);
        }
    }

    private void getMinAndMax() {
        // @TODO: Min and Max optimized
        /*        if (Filename.isStackOrVolume(filename)) {
        ImageDouble img = new ImageDouble();
        img.setFilename(filename);

        double min_max[] = img.getMinAndMax();
        min = min_max[0];
        max = min_max[1];
        } else {*/
        double min_max[] = ImageOperations.getMinAndMax(data);
        min = min_max[0];
        max = min_max[1];
//        }
    }

    public ImagesTableModel(String filenames[], boolean enabled[]) {
        this(filenames);

        if (enabled != null) {
            enableItems(data, enabled);
        }
    }

    public ImagesTableModel(String filenames[]) {
        super();

        loadImages(filenames);
    }

    private String loadStackOrVolume(String filename) {
        try {
            ImageDouble image = new ImageDouble();

            image.readHeader(filename);

            nslices = image.getZsize();
            nimages = image.getNsize();

            for (int i = ImageDouble.FIRST_IMAGE; i <= nimages; i++) {
                for (int j = ImageDouble.FIRST_SLICE; j <= nslices; j++) {
                    addItem(filename, j, i, true);
                }
            }
        } catch (Exception ex) {
            return ex.getMessage();
        }

        return null;
    }

    private void loadImages(String filenames[]) {
        nimages = filenames.length;

        for (int i = 0; i < nimages; i++) {
            String file = Filename.getFilename(filenames[i]);
            long nimage = Filename.getNimage(filenames[i]);

            addItem(file, ImageDouble.FIRST_SLICE, nimage, true);
        }
    }

    private String loadMetaData(String filename) {
        String message = null;

        try {
            md = new MetaData(filename);

            long ids[] = md.findObjects();

            nimages = ids.length;

            for (long id : ids) {
                boolean enabled = true; // If ENABLED field is not present it will be true.
                if (md.containsLabel(MDLabel.MDL_ENABLED)) {
                    enabled = md.getValueInt(MDLabel.MDL_ENABLED, id) == 0 ? false : true;
                }

                String imagefilename = md.getValueString(MDLabel.MDL_IMAGE, id, true);

                long nimage = Filename.getNimage(imagefilename);
                String name = Filename.getFilename(imagefilename);

                addItem(name, ImageDouble.FIRST_SLICE, nimage, enabled);
            }
        } catch (Exception ex) {
            message = ex.getMessage();
            //ex.printStackTrace();
        }

        return message;
    }

    protected void addItem(String filename, int slice, long image, boolean enabled) {
        TableImageItem item = new TableImageItem(new File(filename), slice, image, cache);
        item.setEnabled(enabled);

        data.add(item);
    }

    public String getFilename() {
        return filename;
    }

    public Vector<TableImageItem> getAllItems() {
        return data;
    }

    public LinkedList<TableImageItem> getSelectedItems() {
        return selectedItems;
    }

    @Override
    public String getColumnName(int column) {
        return "" + (column + 1);
    }

    @Override
    public Class getColumnClass(int c) {
        Object object = getValueAt(0, c);

        return object != null ? object.getClass() : null;
    }

    public int getRowCount() {
        return rows;
    }

    public int getColumnCount() {
        return cols;
    }

    public Object getValueAt(int rowIndex, int columnIndex) {
        int index = getDataIndex(rowIndex, columnIndex);

        return index < data.size() ? data.get(index) : null;
    }

    private int getDataIndex(int rowIndex, int columnIndex) {
        return rowIndex * getColumnCount() + columnIndex;
    }

    public int[] getRowColForIndex(int index) {
        int row = index / getColumnCount();
        int col = index % getColumnCount();

        return new int[]{row, col};
    }

    public int getSize() {
        return data.size();
    }

    public String getTitle() {
        String title = "";

        if (filename != null) {
            File f = new File(filename);
            title = f.getName() + " /";
        }

        if (nimages > ImageDouble.FIRST_IMAGE) {
            title += " " + nimages + " images";
        }

        if (nslices > ImageDouble.FIRST_SLICE) {
            title += " " + nslices + " slices";
        }

        return title;
    }

    // Selected items:
    public void setSelected(int row, int col, boolean selected) {
        TableImageItem item = (TableImageItem) getValueAt(row, col);

        if (item != null) {
            item.setSelected(selected);

            if (selected) {
                if (!selectedItems.contains(item)) {
                    selectedItems.addLast(item);
                }
            } else {
                if (selectedItems.contains(item)) {
                    selectedItems.remove(item);
                }
            }
        }
    }

    public void toggleSelected(int row, int col) {
        TableImageItem item = (TableImageItem) getValueAt(row, col);

        if (item != null) {
            setSelected(row, col, !item.isSelected());
        }
    }

    // To select image directly
    public void setSelected(int index) {
        clearSelection();

        TableImageItem item = data.elementAt(index);
        item.setSelected(true);
        selectedItems.addLast(item);
    }

    public void selectAll() {
        clearSelection();

        for (int i = 0; i < getSize(); i++) {
            TableImageItem item = data.elementAt(i);

            item.setSelected(true);
            selectedItems.addLast(item);
        }
    }

    public void invertSelection() {
        for (int i = 0; i < getRowCount(); i++) {
            for (int j = 0; j < getColumnCount(); j++) {
                toggleSelected(i, j);
            }
        }
    }

    public void clearSelection() {
        int n = selectedItems.size();   // Size will change over iterations while index will be increased.

        for (int i = 0; i < n; i++) {
            selectedItems.getLast().setSelected(false);
            selectedItems.removeLast();
        }
    }

    private void enableItems(List<TableImageItem> items, boolean enable[]) {
        if (items != null) {
            for (int i = 0; i < items.size(); i++) {
                items.get(i).setEnabled(enable[i]);
            }
            fireTableDataChanged();
        }
    }

    private void enableItems(List<TableImageItem> items, boolean enable) {
        if (items != null) {
            for (int i = 0; i < items.size(); i++) {
                items.get(i).setEnabled(enable);
            }

            fireTableDataChanged();
        }
    }

    public void enableAllItems() {
        enableItems(data, true);
    }

    public void disableAllItems() {
        enableItems(data, false);
    }

    public void enableSelectedItems() {
        enableItems(selectedItems, true);
    }

    public void disableSelectedItems() {
        enableItems(selectedItems, false);
    }

    public void setZoomScale(double zoomScale) {
        for (int i = 0; i < getSize(); i++) {
            data.elementAt(i).setZoomScale(zoomScale);
        }
    }

    public void setRows(int rows) {
        if (rows > getSize()) {
            rows = getSize();
        }

        if (rows > 0) {
            this.rows = rows;

            cols = getNecessaryCols(rows);

            fireTableStructureChanged();
        }
    }

    public void setColumns(int cols) {
        if (cols > getSize()) {
            cols = getSize();
        }

        if (cols > 0) {
            this.cols = cols;

            rows = getNecessaryRows(cols);

            fireTableStructureChanged();
        }
    }

    public void autoAdjustColumns(int width, int interCellWidth) {
        int displayableColumns = width / (getCellWidth() + 2 * interCellWidth);

        if (getColumnCount() != displayableColumns) {
            setColumns(displayableColumns);
        }
    }

    public int getCellWidth() {
        return getAllItems().elementAt(0).getThumbnailWidth();//(int) (getAllItems().elementAt(0).getWidth() * zoomScale);
    }

    public int getCellHeight() {
        return getAllItems().elementAt(0).getThumbnailHeight();//(int) (getAllItems().elementAt(0).getHeight() * zoomScale);
    }

    @SuppressWarnings("empty-statement")
    private int getNecessaryRows(int cols) {
        int rows_;

        for (rows_ = 1; cols * rows_ < data.size(); rows_++);

        return rows_;
    }

    @SuppressWarnings("empty-statement")
    private int getNecessaryCols(int rows) {
        int cols_;

        for (cols_ = 1; cols_ * rows < data.size(); cols_++);

        return cols_;
    }

    public void setNormalized(boolean normalize) {
        this.normalize = normalize;

        if (min == Double.POSITIVE_INFINITY && max == Double.NEGATIVE_INFINITY) {
            DEBUG.printMessage(" +++ Retrieving min and max.");

            getMinAndMax();
        }

        DEBUG.printMessage(" >>> Normalize " + (normalize ? "ON" : "OFF") + " > "
                + (normalize ? "m=" + min + "/M=" + max : ""));
    }

    public boolean isNormalizing() {
        return min > Double.NEGATIVE_INFINITY && max < Double.POSITIVE_INFINITY && normalize;
    }

    public double getNormalizeMin() {
        return min;
    }

    public double getNormalizeMax() {
        return max;
    }

    public boolean isStack() {
        return filename != null && Filename.isStack(filename);
    }

    public boolean isVolume() {
        return filename != null && Filename.isVolume(filename);
    }
}
