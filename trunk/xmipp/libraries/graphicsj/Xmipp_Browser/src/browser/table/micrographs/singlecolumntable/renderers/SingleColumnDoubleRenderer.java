/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs.singlecolumntable.renderers;

import browser.table.micrographs.renderers.MicrographDoubleRenderer;
import javax.swing.JTable;

/**
 *
 * @author Juanjo Vega
 */
public class SingleColumnDoubleRenderer extends MicrographDoubleRenderer {

    @Override
    protected boolean isRowEnabled(JTable table, int row) {
        return true;
    }
}
