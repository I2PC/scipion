/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package metadata.renderers;

import java.awt.Component;
import java.awt.Font;
import javax.swing.JLabel;
import javax.swing.JTable;
import metadata.images.TableFileItem;

/**
 *
 * @author Juanjo Vega
 */
public class FileItemRenderer extends MetaDataRowDisablerRenderer {

    private Font font;

    public FileItemRenderer() {
        super();

        setHorizontalAlignment(JLabel.CENTER);

        font = getFont();
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        TableFileItem item = (TableFileItem) value;

        String str = getShortLabel(item.getOriginalValue(), table.getColumnModel().getColumn(column).getWidth());

        setToolTipText(item.getPath());

        return super.getTableCellRendererComponent(table, str, isSelected, hasFocus, row, column);
    }

    private String getShortLabel(String label, int width) {
        StringBuffer sb = new StringBuffer(label);
        String sortLabel = sb.toString();

        int w = getFontMetrics(font).stringWidth(sortLabel);

        int i = 0;
        while (w > width) {
            sortLabel = "..." + sb.substring(i++);
            w = getFontMetrics(font).stringWidth(sortLabel);
        }

        return sortLabel;
    }
}
