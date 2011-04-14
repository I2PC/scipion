/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs.renderers;

import browser.imageitems.TableImageItem;
import java.awt.Component;
import java.awt.Font;
import javax.swing.JLabel;
import javax.swing.JTable;

/**
 *
 * @author Juanjo Vega
 */
public class MicrographFileNameRenderer extends MicrographRowDisablerRenderer {

    public MicrographFileNameRenderer() {
        super();

        setHorizontalAlignment(JLabel.CENTER);
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        TableImageItem item = (TableImageItem) value;

        StringBuffer sb = new StringBuffer(item.getFileName());
        String str = sb.toString();

        int columnWidth = table.getColumnModel().getColumn(column).getWidth();

        Font font = getFont();
        int w = getFontMetrics(font).stringWidth(str);

        int i = 0;
        while (w > columnWidth) {
            str = "..." + sb.substring(i++);
            w = getFontMetrics(font).stringWidth(str);
        }

        setToolTipText(item.getFileName());

        return super.getTableCellRendererComponent(table, str, isSelected, hasFocus, row, column);//this;
    }
}
