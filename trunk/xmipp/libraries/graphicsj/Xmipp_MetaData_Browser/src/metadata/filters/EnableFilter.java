/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package metadata.filters;

import metadata.models.MetaDataTableModel;
import javax.swing.RowFilter;
import javax.swing.RowFilter.Entry;

/**
 *
 * @author Juanjo Vega
 */
public class EnableFilter extends RowFilter<MetaDataTableModel, Object> {

    protected boolean filtering = false;

    public void setFiltering(boolean filtering) {
        this.filtering = filtering;
    }

    public boolean isFiltering() {
        return filtering;
    }

    @Override
    public boolean include(Entry<? extends MetaDataTableModel, ? extends Object> entry) {
        return !filtering || entry.getValue(entry.getModel().getEnabledColumnIndex()).equals(true);
    }
}
