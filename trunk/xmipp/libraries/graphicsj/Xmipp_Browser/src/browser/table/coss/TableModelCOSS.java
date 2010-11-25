/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.coss;

import java.awt.Image;
import java.io.File;
import javax.swing.RowSorter;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;

/**
 *
 * @author Juanjo Vega
 */
public class TableModelCOSS extends DefaultTableModel {

    private String columnsNames[] = new String[]{
        "Selected", "Image 1", "Image 2", "Image 3", "Info"
    };
    private Class columnsClasses[] = {
        Boolean.class,
        Image.class,
        Image.class,
        Image.class,
        String.class};
    private RowSorter rowSorter;
    private JTableHeader header;

    public TableModelCOSS() {
        super();

        for (int i = 0; i < columnsNames.length; i++) {
            addColumn(columnsNames[i]);
        }
    }

    public void setRowSorter(RowSorter rowSorter) {
        this.rowSorter = rowSorter;
    }

    public void setHeader(JTableHeader header) {
        this.header = header;
    }

    public void appendEmptyRow() {
        Object[] newRow = new Object[]{false, "", "", "", ""};

        int row = getRowCount();

        addRow(newRow);

        fireTableRowsInserted(row, row);
    }

    /*
     * JTable uses this method to determine the default renderer/
     * editor for each cell.  If we didn't implement this method,
     * then the first column would contain text ("true"/"false"),
     * rather than a check box.
     */
    @Override
    public Class getColumnClass(int c) {
        return columnsClasses[c];
    }

    @Override
    public void setValueAt(Object aValue, int row, int column) {
        super.setValueAt(aValue, row, column);

        fireTableRowsUpdated(row, row); // Update the entire row, so related cells will show changes.
    }

    public void save(File file) {
        System.out.println(" ### Saving: " + file);

        for (int col = 0; col < getColumnCount(); col++) {
            System.out.print(header.getColumnModel().getColumn(col).getModelIndex() + "\t");
        }
        System.out.println();

        for (int row = 0; row < getRowCount(); row++) {
            int row_ = rowSorter.convertRowIndexToModel(row);
            for (int col = 0; col < getColumnCount(); col++) {
                int col_ = header.getColumnModel().getColumn(col).getModelIndex();
                System.out.print(getValueAt(row_, col_) + "\t");
            }
            System.out.println();
        }
    }
}
