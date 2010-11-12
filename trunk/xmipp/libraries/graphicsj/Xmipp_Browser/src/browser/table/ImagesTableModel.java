/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table;

import browser.Cache;
import browser.files.FilterFilesModel;
import browser.imageitems.TableImageItem;
import browser.imageitems.listitems.AbstractImageItem;
import ij.ImagePlus;
import java.io.File;
import java.util.Vector;
import javax.swing.table.AbstractTableModel;
import xmipp.io.Sel_Reader;
import xmipp.io.Spider_Reader;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesTableModel extends AbstractTableModel {

    private Vector<TableImageItem> data = new Vector<TableImageItem>();
    private Vector<TableImageItem> selectedItems = new Vector<TableImageItem>();
    private int rows, cols;
    private boolean showLabels = true;
    private Cache cache = new Cache();
    protected boolean normalize;
    protected double min, max;

    public ImagesTableModel() {
        this(1, 1);
    }

    public ImagesTableModel(int rows, int cols) {
        super();

        this.rows = rows;
        this.cols = cols;
    }

    /**
     * Refreshes entire table by clearing cache.
     */
    public void refresh() {
        cache.clear();

//        updateNormalizeValues();
    }

    /**
     * Clears data.
     */
    public void clear() {
        clearSelection();

        data.clear();
        rows = cols = 1;
    }

    public void setNormalized(double min, double max) {
        this.min = min;
        this.max = max;

        setNormalized(true);
    }

    public void setNormalizedAuto() {
        double min_max[] = ImageTableOperations.getMinAndMax(data);

        setNormalized(min_max[0], min_max[1]);
    }

    public void disableNormalization() {
        setNormalized(false);
    }

    private void setNormalized(boolean normalize) {
        this.normalize = normalize;

        // Refreshes items.
        for (int i = 0; i < data.size(); i++) {
            data.get(i).setNormalized(normalize);
        }
    }

    public void toggleSelected(int row, int col) {
        TableImageItem item = (TableImageItem) getValueAt(row, col);

        if (item != null) {
            setSelected(row, col, !item.isSelected());
        }
    }

    public void setSelected(int row, int col, boolean selected) {
        TableImageItem item = (TableImageItem) getValueAt(row, col);

        if (item != null) {
            item.setSelected(selected);

            if (selected) {
                if (!selectedItems.contains(item)) {
                    selectedItems.add(item);
                }
            } else {
                if (selectedItems.contains(item)) {
                    selectedItems.remove(item);
                }
            }
        }
    }

    public Vector<TableImageItem> getItems() {
        return data;
    }

    public Vector<TableImageItem> getSelectedItems() {
        return selectedItems;
    }

    public Vector<TableImageItem> getSelectedAndEnabledItems() {
        Vector<TableImageItem> enabledItems = new Vector<TableImageItem>();

        for (int i = 0; i < selectedItems.size(); i++) {
            TableImageItem tableImageItem = selectedItems.elementAt(i);
            if (tableImageItem.isEnabled()) {
                enabledItems.add(tableImageItem);
            }
        }

        return enabledItems;
    }

    public double getMin() {
        //System.out.println("m=" + min);
        return min;
    }

    public double getMax() {
        //System.out.println("M=" + max);
        return max;
    }

    public void setSelected(int index) {
        clearSelection();

        TableImageItem item = data.elementAt(index);
        item.setSelected(true);
        selectedItems.add(item);
    }

    public void selectAll() {
        clearSelection();

        for (int i = 0; i < getSize(); i++) {
            TableImageItem item = data.elementAt(i);

            item.setSelected(true);
            selectedItems.add(item);
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
        for (int i = 0; i < selectedItems.size(); i++) {
            selectedItems.elementAt(i).setSelected(false);
        }

        selectedItems.clear();
    }

    public void addImage(AbstractImageItem itemImage) {
        data.add(new TableImageItem(itemImage, cache, this));

        fireTableStructureChanged();
    }

    public void addVolume(AbstractImageItem itemImage) {
        for (int i = 0; i < itemImage.nslices; i++) {
            data.add(new TableImageItem(itemImage, cache, i + 1, this));
        }

        fireTableStructureChanged();
    }

    public int addVolumeFile(File file) {
        String name = file.getName();
        String parent = file.getParent();
        if (parent == null) {
            parent = System.getProperty("user.dir");
        }

        ImagePlus ip = (new Spider_Reader()).load(parent, name);

        File fullFile = new File(parent + File.separator + name);
        for (int i = 0; i < ip.getStackSize(); i++) {
            data.add(new TableImageItem(
                    (AbstractImageItem) FilterFilesModel.createSuitableFileItem(fullFile),
                    cache, i + 1, this));
        }

        fireTableStructureChanged();

        return ip.getStackSize();
    }

    public void addSelVolume(String dir, String file) {
        String fileNames[] = Sel_Reader.loadFileNames(dir, file);

        for (int i = 0; i < fileNames.length; i++) {
            data.add(new TableImageItem(
                    (AbstractImageItem) FilterFilesModel.createSuitableFileItem(new File(fileNames[i])),
                    cache, this));
        }

        fireTableStructureChanged();
    }

    @Override
    public String getColumnName(int column) {
        return "" + (column + 1);
    }

    public int getRowCount() {
        return rows;
    }

    public int getColumnCount() {
        return cols;
    }

    public void setRows(int rows) {
        if (rows > 0 && rows <= getSize()) {
            this.rows = rows;

            cols = getNecessaryCols(rows);

            fireTableStructureChanged();
        }
    }

    public void setColumns(int cols) {
        if (cols > 0 && cols <= getSize()) {
            this.cols = cols;

            rows = getNecessaryRows(cols);

            fireTableStructureChanged();
        }
    }

    public void setShowLabels(boolean show) {
        showLabels = show;
    }

    public boolean isShowingLabels() {
        return showLabels;
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

    /**
     * Returns index at data vector depending on the specified items arrangement.
     * @param rowIndex
     * @param columnIndex
     * @return
     */
    public int getDataIndex(int rowIndex, int columnIndex) {
        return rowIndex * getColumnCount() + columnIndex;
    }

    public int[] getRowColForIndex(int index) {
        int row = index / getColumnCount();
        int col = index % getColumnCount();

        return new int[]{row, col};
    }

    public Object getValueAt(int rowIndex, int columnIndex) {
        int index = getDataIndex(rowIndex, columnIndex);

        return index < data.size() ? data.get(index) : null;
    }

    public int getSize() {
        return data.size();
    }

    @Override
    public Class getColumnClass(int c) {
        Object object = getValueAt(0, c);

        return object != null ? object.getClass() : null;
    }

    /**
     * Sets zoom in all objects.
     * @param zoom
     */
    public void setZoom(int zoom) {
        for (int i = 0; i < data.size(); i++) {
            data.elementAt(i).setZoom(zoom);
        }
    }

    protected void resizeCache(int newLimit) {
        cache.resize(newLimit);
    }
}
