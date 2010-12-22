/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs.renderers;

import browser.table.micrographs.TableImageItemCOSS;
import java.awt.Component;
import java.awt.Font;
import javax.swing.JLabel;
import javax.swing.JTable;

/**
 *
 * @author Juanjo Vega
 */
public class FileNameRenderer extends DisablerRenderer {

    public FileNameRenderer() {
        super();

        label.setHorizontalAlignment(JLabel.CENTER);
    }

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        TableImageItemCOSS item = (TableImageItemCOSS) value;

        StringBuffer sb = new StringBuffer(item.getLabel());
        String s = item.getLabel();

        int columnWidth = table.getColumnModel().getColumn(column).getWidth();

        Font font = label.getFont();
        int w = label.getFontMetrics(font).stringWidth(s);

        int i = 0;
        while (w > columnWidth) {
            s = "..." + sb.substring(i++);
            w = label.getFontMetrics(font).stringWidth(s);
        }

        label.setText("<html><font color=\"" + getForegroundColor(table, row) + "\">" + s + "</font></html>");
        setToolTipText(item.getFileName());

        setCellBackgroundForeground(table, isSelected, hasFocus);

        return this;
    }
}
