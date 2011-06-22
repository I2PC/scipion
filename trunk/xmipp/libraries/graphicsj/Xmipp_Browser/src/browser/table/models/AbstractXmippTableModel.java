/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.models;

import browser.Cache;
import browser.imageitems.tableitems.AbstractTableImageItem;
import ij.IJ;
import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;
import javax.swing.table.AbstractTableModel;

/**
 *
 * @author Juanjo Vega
 */
public abstract class AbstractXmippTableModel extends AbstractTableModel {

    protected String filename;
    protected Vector<AbstractTableImageItem> data = new Vector<AbstractTableImageItem>();
    protected LinkedList<AbstractTableImageItem> selectedItems = new LinkedList<AbstractTableImageItem>();
    protected int rows, cols;
    protected Cache cache = new Cache();
    protected boolean normalize = false;
    protected double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;

    public AbstractXmippTableModel() {
        super();
    }

    public AbstractXmippTableModel(String filename) {
        super();

        this.filename = filename;

        String message = null;

        File f = new File(filename);
        if (f.exists()) {
            message = populateTable(filename);
        } else {
            message = "File not found: " + filename;
        }

        if (message != null) {
            IJ.error(message);
        }

//        printMD(md);
    }

    protected abstract String populateTable(String filename);

//    private static void printMD(MetaData md) {
//        int labels[] = md.getActiveLabels();
//
//        for (int i : labels) {
//            String label = MetaData.label2Str(i);
//            System.out.println(i + ": " + label);
//        }
//
//        long ids[] = md.findObjects();
//        for (long id : ids) {
//            System.out.println(id + ": "
//                    + md.getValueString(MDLabel.MDL_IMAGE, id, true) + " / "
//                    + md.getValueInt(MDLabel.MDL_ENABLED, id));
//        }
//
//        System.out.println("++++++++++++++++++++++++++++++++++++++++");
//    }
    protected abstract void getMinAndMax();

    public abstract String getFilename();

    public Vector<AbstractTableImageItem> getAllItems() {
        return data;
    }

    public LinkedList<AbstractTableImageItem> getSelectedItems() {
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

    public abstract String getTitle();

    // Selected items:
    public void setSelected(int row, int col, boolean selected) {
        AbstractTableImageItem item = (AbstractTableImageItem) getValueAt(row, col);

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
        AbstractTableImageItem item = (AbstractTableImageItem) getValueAt(row, col);

        if (item != null) {
            setSelected(row, col, !item.isSelected());
        }
    }

    // To select image directly
    public void setSelected(int index) {
        clearSelection();

        AbstractTableImageItem item = data.elementAt(index);
        item.setSelected(true);
        selectedItems.addLast(item);
    }

    public void selectAll() {
        clearSelection();

        for (int i = 0; i < getSize(); i++) {
            AbstractTableImageItem item = data.elementAt(i);

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
//
//    private void enableItems(List<AbstractTableImageItem> items, boolean enable[]) {
//        if (items != null) {
//            for (int i = 0; i < items.size(); i++) {
//                items.get(i).setEnabled(enable[i]);
//            }
//            fireTableDataChanged();
//        }
//    }

    private void enableItems(List<AbstractTableImageItem> items, boolean enable) {
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

            // Calculates necessary cols for desired number of rows.
            cols = (int) Math.ceil((double) data.size() / (double) rows);

//            DEBUG.printMessage("R:" + rows + " / C:" + cols + " > S: " + (rows * cols));

            fireTableStructureChanged();
        }
    }

    public void setColumns(int cols) {
        if (cols > getSize()) {
            cols = getSize();
        }

        if (cols > 0) {
            this.cols = cols;

            // Calculates necessary rows for desired number of cols.
            rows = (int) Math.ceil((double) data.size() / (double) cols);

//            DEBUG.printMessage("R:" + rows + " / C:" + cols + " > S: " + (rows * cols));

            fireTableStructureChanged();
        }
    }

    public void autoAdjustColumns(int width, int interCellWidth) {
        int displayableColumns = (int) Math.floor(
                width / (getCellWidth() + 2 * interCellWidth));

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

    public void setNormalized(boolean normalize) {
        this.normalize = normalize;

        if (min == Double.POSITIVE_INFINITY && max == Double.NEGATIVE_INFINITY) {
            getMinAndMax();
        }

        /*        DEBUG.printMessage(" >>> Normalize " + (normalize ? "ON" : "OFF") + " > "
        + (normalize ? "m=" + min + "/M=" + max : ""));*/
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

    public abstract boolean isStack();

    public abstract boolean isVolume();
}
