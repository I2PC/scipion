/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.metadata.models;

import javax.swing.JTable;
import javax.swing.ListModel;
import javax.swing.event.ListDataListener;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesRowHeaderModel implements ListModel {

    private JTable table;

    public ImagesRowHeaderModel(JTable table) {
        super();

        this.table = table;
    }

    public int getSize() {
        return table.getRowCount();
    }

    public Object getElementAt(int i) {
        return i;
    }

    public void addListDataListener(ListDataListener ll) {
    }

    public void removeListDataListener(ListDataListener ll) {
    }
}
