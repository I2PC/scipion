/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.renderers;

import java.awt.Component;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.table.TableCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public abstract class RowDisablerRenderer extends JLabel implements TableCellRenderer {

    protected Border FocusedBorder = UIManager.getBorder("Table.focusCellHighlightBorder");
    private int ENABLED_COLUMN_INDEX = 0;

    public RowDisablerRenderer() {
        super();

        setHorizontalAlignment(JLabel.LEFT);
        setVerticalAlignment(JLabel.CENTER);
    }

    abstract public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column);

    protected boolean isRowEnabled(JTable table, int row) {
        return (Boolean) table.getValueAt(row, ENABLED_COLUMN_INDEX);
    }

    protected String getForegroundColor(JTable table, int row) {
        String color = isRowEnabled(table, row) ? "#0A0A0A" : "#A9A9A9";

        return color;
    }

    protected void setCellBackgroundForeground(JTable table, boolean isSelected, boolean hasFocus) {
        if (isSelected) {
            setBackground(table.getSelectionBackground());
            setForeground(table.getSelectionForeground());
        } else {
            setBackground(table.getBackground());
            setForeground(table.getForeground());
        }

        if (hasFocus) {
            setBorder(FocusedBorder);
        } else {
            setBorder(null);
        }
    }
}
