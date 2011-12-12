/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package metadata.renderers;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;
import metadata.models.XmippTableModelRowDisabler;

/**
 *
 * @author Juanjo Vega
 */
public class MetaDataRowDisablerRenderer extends DefaultTableCellRenderer {

    private Font font;

    public MetaDataRowDisablerRenderer() {
        super();

        setHorizontalAlignment(JLabel.LEFT);
        setVerticalAlignment(JLabel.CENTER);


        font = getFont();
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);

        if (value != null) {
//            GalleryImageItem item = (GalleryImageItem) value;

            String labelStr = value.toString();//item.getOriginalValue();
            String columnName = table.getColumnName(column);

            int length = Math.max(labelStr.length(), columnName.length());

            int w = font.getSize() * length;
            int h = font.getSize() * 3 / 2;

            String label = getShortLabel(/*item.getOriginalValue()*/labelStr, w);

            setText("<html><font color=\"" + getFontForegroundColor(table, row) + "\">" + label + "</font></html>");

            setPreferredSize(new Dimension(w, h));
        }

        return this;
    }

    protected String getShortLabel(String label, int width) {
        StringBuilder sb = new StringBuilder(label);
        String sortLabel = sb.toString();

        int w = getFontMetrics(font).stringWidth(sortLabel);

        int i = 0;
        while (w > width) {
            sortLabel = "..." + sb.substring(i++);
            w = getFontMetrics(font).stringWidth(sortLabel);
        }

        return sortLabel;
    }

    protected boolean isRowEnabled(JTable table, int row) {
        int modelRow = table.convertRowIndexToModel(row);

        return ((XmippTableModelRowDisabler) table.getModel()).isRowEnabled(modelRow);
    }

    protected String getFontForegroundColor(JTable table, int row) {
        String color = isRowEnabled(table, row) ? "#0A0A0A" : "#A9A9A9";

        return color;
    }
}
