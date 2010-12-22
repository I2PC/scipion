/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs.renderers;

import java.awt.BorderLayout;
import java.awt.Component;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.table.TableCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public abstract class DisablerRenderer extends JPanel implements TableCellRenderer {

    JLabel label = new JLabel();
    Border FocusedBorder = UIManager.getBorder("Table.focusCellHighlightBorder");

    public DisablerRenderer() {
        super();

        setLayout(new BorderLayout());
        add(label, BorderLayout.CENTER);

        label.setHorizontalAlignment(JLabel.LEFT);
        label.setVerticalAlignment(JLabel.CENTER);
    }

    abstract public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column);

    protected boolean isRowEnabled(JTable table, int row) {
        return (Boolean) table.getValueAt(row, 0);
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
