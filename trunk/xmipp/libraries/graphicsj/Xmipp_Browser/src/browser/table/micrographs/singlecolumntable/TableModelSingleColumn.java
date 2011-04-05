/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs.singlecolumntable;

import browser.table.micrographs.MicrographsTableModel;
import java.util.Arrays;
import javax.swing.table.AbstractTableModel;

/**
 *
 * @author Juanjo Vega
 */
public class TableModelSingleColumn extends AbstractTableModel {

    private MicrographsTableModel externTableModel;
    private int column;

    public TableModelSingleColumn(MicrographsTableModel externTableModel, int column) {
        this.externTableModel = externTableModel;
        this.column = column;
    }

    public int getRowCount() {
        return externTableModel.getRowCount();
    }

    public int getColumnCount() {
        return 1;
    }

    public Object getValueAt(int row, int column) {
        return externTableModel.getValueAt(row, this.column);
    }

    @Override
    public String getColumnName(int column) {
        return externTableModel.getColumnName(this.column);
    }

    @Override
    public Class<?> getColumnClass(int columnIndex) {
        return externTableModel.getColumnClass(this.column);
    }

    public boolean isImage() {
        return Arrays.binarySearch(MicrographsTableModel.imagesColumnIndex, column) >= 0;
    }

    public boolean isFilename() {
        return Arrays.binarySearch(MicrographsTableModel.filenameColumnIndex, column) >= 0;
    }

    public boolean isDouble() {
        return Arrays.binarySearch(MicrographsTableModel.doubleColumnIndex, column) >= 0;
    }
}
