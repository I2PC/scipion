/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package micrographs.renderers;

import browser.table.micrographs.MicrographsTableModel;
import java.awt.Component;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public class MicrographRowDisablerRenderer extends DefaultTableCellRenderer {

    public MicrographRowDisablerRenderer() {
        super();

        setHorizontalAlignment(JLabel.LEFT);
        setVerticalAlignment(JLabel.CENTER);
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        String str = "";
        if (value != null) {
            str = "<html><font color=\"" + getFontForegroundColor(table, row) + "\">" + value + "</font></html>";
        }

        return super.getTableCellRendererComponent(table, str, isSelected, hasFocus, row, column);
    }

    protected boolean isRowEnabled(JTable table, int row) {
        int modelRow = table.convertRowIndexToModel(row);

        return ((MicrographsTableModel) table.getModel()).isRowEnabled(modelRow);
    }

    protected String getFontForegroundColor(JTable table, int row) {
        String color = isRowEnabled(table, row) ? "#0A0A0A" : "#A9A9A9";

        return color;
    }
}
