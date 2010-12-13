/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.coss.renderers;

import browser.table.coss.TableImageItemCOSS;
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
public class FileNameRenderer extends JPanel implements TableCellRenderer {

    JLabel label = new JLabel();

    public FileNameRenderer() {
        super();

        setLayout(new BorderLayout());
        add(label, BorderLayout.CENTER);

        label.setHorizontalAlignment(JLabel.CENTER);
        label.setVerticalAlignment(JLabel.CENTER);
    }

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        TableImageItemCOSS item = (TableImageItemCOSS) value;

        StringBuffer sb = new StringBuffer(item.getLabel());
        String s;

        int length = sb.length();
        if (length > 10) {
            s = "..." + sb.substring(length - 10, length);
        } else {
            s = sb.toString();
        }

        label.setText(s);
        setToolTipText(item.getFileName());

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
