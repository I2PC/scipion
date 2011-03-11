/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.renderers;

import browser.imageitems.TableImageItem;
import java.awt.Component;
import java.awt.Font;
import javax.swing.JLabel;
import javax.swing.JTable;

/**
 *
 * @author Juanjo Vega
 */
public class FileNameRenderer extends RowDisablerRenderer {

    public FileNameRenderer() {
        super();

        setHorizontalAlignment(JLabel.CENTER);
    }

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        TableImageItem item = (TableImageItem) value;

        StringBuffer sb = new StringBuffer(item.getLabel());
        String s = item.getLabel();

        int columnWidth = table.getColumnModel().getColumn(column).getWidth();

        Font font = getFont();
        int w = getFontMetrics(font).stringWidth(s);

        int i = 0;
        while (w > columnWidth) {
            s = "..." + sb.substring(i++);
            w = getFontMetrics(font).stringWidth(s);
        }

        setText("<html><font color=\"" + getForegroundColor(table, row) + "\">" + s + "</font></html>");
        setToolTipText(item.getFileName());

        setCellBackgroundForeground(table, isSelected, hasFocus);

        return this;
    }
}
