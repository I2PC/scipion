/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.renderers;

import java.awt.Component;
import javax.swing.JLabel;
import javax.swing.JTable;

/**
 *
 * @author Juanjo Vega
 */
public class DoubleRenderer extends RowDisablerRenderer {

    public DoubleRenderer() {
        super();

        setHorizontalAlignment(JLabel.RIGHT);
    }

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        Double double_ = (Double) value;

        setText("<html><font color=\"" + getForegroundColor(table, row) + "\">" + double_ + "</font></html>");

        setCellBackgroundForeground(table, isSelected, hasFocus);

        return this;
    }
}
