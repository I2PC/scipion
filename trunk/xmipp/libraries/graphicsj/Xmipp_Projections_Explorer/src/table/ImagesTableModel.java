/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package table;

import ij.ImagePlus;
import java.util.Vector;
import javax.swing.table.AbstractTableModel;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesTableModel extends AbstractTableModel {

    protected Vector<ScoreItem> data = new Vector<ScoreItem>();
    protected ImagePlus focusedItem;

    public ImagesTableModel() {
        super();
    }

    public void clear() {
        data.clear();
    }

    public void addItem(ScoreItem scoreItem) {
        data.add(scoreItem);

        fireTableStructureChanged();
    }

    @Override
    public String getColumnName(int column) {
        return "" + column;
    }

    public int getRowCount() {
        return 1;
    }

    public int getColumnCount() {
        return data.size();
    }

    public Object getValueAt(int rowIndex, int columnIndex) {
        return columnIndex < data.size() ? data.get(columnIndex) : null;
    }

    @Override
    public Class getColumnClass(int c) {
        Object object = getValueAt(0, c);

        return object != null ? object.getClass() : null;
    }

    public Vector<ScoreItem> getData() {
        return data;
    }
}
