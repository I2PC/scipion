/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.rotspectra;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.LinkedList;
import javax.swing.table.AbstractTableModel;
import xmipp.jni.Filename;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class RotSpectraTableModel extends AbstractTableModel {

    final static String blockKerDenSOM = "KerDenSOM_Layout" + Filename.SEPARATOR;
    final static String blockVectorContent = "vectorContent" + Filename.SEPARATOR;
    int w, h;
    int vectorsSize;
    LinkedList<RotSpectraVector> data = new LinkedList<RotSpectraVector>();
    protected LinkedList<RotSpectraVector> selectedItems = new LinkedList<RotSpectraVector>();
    String filenameClasses, filenameVectors, filenameData;
    boolean showLabels;

    public RotSpectraTableModel(String filenameClasses, String filenameVectors, String filenameData) throws Exception {
        this.filenameClasses = filenameClasses;
        this.filenameVectors = filenameVectors;
        this.filenameData = filenameData;

        load(filenameClasses, filenameVectors, filenameData);
    }

    void load(String filenameClasses, String filenameVectors, String filenameData) throws Exception {
        MetaData mdClasses = new MetaData(blockKerDenSOM + filenameClasses);
        MetaData mdVectors = new MetaData(filenameVectors);

        long id = mdClasses.firstObject();
        w = mdClasses.getValueInt(MDLabel.MDL_XSIZE, id);
        h = mdClasses.getValueInt(MDLabel.MDL_YSIZE, id);
        int nvectors = w * h;

        String blocks[] = getBlocks(blockVectorContent + filenameVectors);

        id = mdVectors.firstObject();
        vectorsSize = mdVectors.getValueInt(MDLabel.MDL_CLASSIFICATION_DATA_SIZE, id);

        double vectors[][] = loadVectors(filenameData, nvectors, vectorsSize);

        for (int i = 0; i < vectors.length; i++) {
            data.add(new RotSpectraVector(
                    blocks[i], filenameClasses, vectors[i]));
        }
    }

    static String[] getBlocks(String filename) throws Exception {
        MetaData mdVectors = new MetaData(filename);

        long ids[] = mdVectors.findObjects();
        String blocks[] = new String[ids.length];
        for (int i = 0; i < ids.length; i++) {
            blocks[i] = mdVectors.getValueString(MDLabel.MDL_IMAGE, ids[i]);
        }

        return blocks;
    }

    static double[][] loadVectors(String filename, int nvectors, int size) throws Exception {
        double vectors[][] = new double[nvectors][size];

        // Load vectors.
        DataInputStream in = new DataInputStream(
                new FileInputStream(filename));

        for (int i = 0; i < nvectors; i++) {
            for (int j = 0; j < size; j++) {
                vectors[i][j] = in.readFloat();
            }
        }

        in.close();

        return vectors;
    }

    @Override
    public String getColumnName(int i) {
        return String.valueOf(i + 1);
    }

    public int getSize() {
        return data.size();
    }

    public int getRowCount() {
        return w;
    }

    public int getColumnCount() {
        return h;
    }

    @Override
    public Class getColumnClass(int i) {
        return getSize() > 0 ? getValueAt(0, 0).getClass() : Object.class;
    }

    public Object getValueAt(int row, int column) {
        int index = row * getColumnCount() + column;

        return data.get(index);
    }

    public void setShowLabels(boolean showLabels) {
        this.showLabels = showLabels;
    }

    public boolean isShowingLabels() {
        return showLabels;
    }

    public void selectRange(int initial_row, int initial_col, int final_row, int final_col) {
        int row0 = Math.min(initial_row, final_row);
        int row1 = Math.max(initial_row, final_row);
        int col0 = Math.min(initial_col, final_col);
        int col1 = Math.max(initial_col, final_col);

        for (int i = row0; i <= row1; i++) {
            for (int j = col0; j <= col1; j++) {
                setSelected(i, j, true);
            }
        }
    }

    public void setSelected(int row, int col, boolean selected) {
        RotSpectraVector item = (RotSpectraVector) getValueAt(row, col);

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
        //fireTableCellUpdated(row, col);
    }

    public ArrayList<RotSpectraVector> getSelectedItems() {
        return new ArrayList<RotSpectraVector>(selectedItems);
    }

    public void toggleSelected(int row, int col) {
        RotSpectraVector item = (RotSpectraVector) getValueAt(row, col);

        if (item != null) {
            setSelected(row, col, !item.isSelected());
        }
    }

    // To select image directly
    public void setSelected(int index) {
        clearSelection();

        RotSpectraVector item = data.get(index);
        item.setSelected(true);
        selectedItems.addLast(item);
    }

    public void selectAll() {
        clearSelection();

        for (int i = 0; i < getSize(); i++) {
            RotSpectraVector item = data.get(i);

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
}
