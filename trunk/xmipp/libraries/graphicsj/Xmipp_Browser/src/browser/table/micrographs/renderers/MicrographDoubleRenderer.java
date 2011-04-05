/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs.renderers;

import browser.table.micrographs.renderers.MicrographRowDisablerRenderer;
import java.awt.Component;
import javax.swing.JLabel;
import javax.swing.JTable;

/**
 *
 * @author Juanjo Vega
 */
public class MicrographDoubleRenderer extends MicrographRowDisablerRenderer {

    public MicrographDoubleRenderer() {
        super();

        setHorizontalAlignment(JLabel.RIGHT);
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        return super.getTableCellRendererComponent(table, "<b>" + value + "</b>", isSelected, hasFocus, row, column);//this;
    }
}
