/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.metadata.filters;

import xmipp.viewer.metadata.models.MetaDataTableModel;
import javax.swing.RowFilter;
import javax.swing.RowFilter.Entry;

/**
 *
 * @author Juanjo Vega
 */
public class RowEnableFilter extends RowFilter<MetaDataTableModel, Object> {

    @Override
    public boolean include(Entry<? extends MetaDataTableModel, ? extends Object> entry) {
        return entry.getValue(entry.getModel().getEnabledColumnIndex()).equals(true);
    }
}
