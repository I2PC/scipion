/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table;

import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.TableColumn;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesTableColumnModel extends DefaultTableColumnModel {

    protected int width;

    @Override
    public void addColumn(TableColumn tc) {
        tc.setPreferredWidth(width);
        super.addColumn(tc);
    }

    public void setWidth(int width) {
        this.width = width;

        autoadjustColumns();
    }

    public void autoadjustColumns() {
        for (int i = 0; i < getColumnCount(); i++) {
            getColumn(i).setPreferredWidth(width);
        }
    }
}
