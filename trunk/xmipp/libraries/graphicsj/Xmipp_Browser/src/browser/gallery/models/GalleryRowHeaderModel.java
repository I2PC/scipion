/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.gallery.models;

import javax.swing.JTable;
import javax.swing.ListModel;
import javax.swing.event.ListDataListener;

/**
 *
 * @author Juanjo Vega
 */
public class GalleryRowHeaderModel implements ListModel {

    private JTable table;
    private int first_index = 0;

    public GalleryRowHeaderModel(JTable table, int first_index) {
        this(table);
        this.first_index = first_index;
    }

    public GalleryRowHeaderModel(JTable table) {
        super();

        this.table = table;
    }

    public int getSize() {
        return table.getRowCount();
    }

    public Object getElementAt(int i) {
        return new Integer(first_index + i);
    }

    public void addListDataListener(ListDataListener ll) {
    }

    public void removeListDataListener(ListDataListener ll) {
    }
}
