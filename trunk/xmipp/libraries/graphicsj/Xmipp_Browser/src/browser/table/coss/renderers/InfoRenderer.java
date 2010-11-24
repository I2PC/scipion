/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.coss.renderers;

import java.awt.BorderLayout;
import java.awt.Component;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.UIManager;
import javax.swing.table.TableCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public class InfoRenderer extends JPanel implements TableCellRenderer {

    JLabel label = new JLabel();

    public InfoRenderer() {
        super();

        setLayout(new BorderLayout());
        add(label, BorderLayout.CENTER);

        label.setHorizontalAlignment(JLabel.LEFT);
        label.setVerticalAlignment(JLabel.CENTER);
    }

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        String string = (String) value;

        boolean selected = (Boolean) table.getModel().getValueAt(row, 0);

        String color = selected ? "#0A0A0A" : "#A9A9A9";

        label.setText("<html><font color=\"" + color + "\">" + string + "</font></html>");

        if (isSelected) {
            setBackground(table.getSelectionBackground());
            setForeground(table.getSelectionForeground());
        } else {
            setBackground(table.getBackground());
            setForeground(table.getForeground());
        }

        if (hasFocus) {
            setBorder(UIManager.getBorder("Table.focusCellHighlightBorder"));
        } else {
            setBorder(null);
        }

        return this;
    }
}
