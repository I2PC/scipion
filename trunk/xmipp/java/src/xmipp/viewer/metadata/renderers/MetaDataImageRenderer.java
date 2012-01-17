/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.metadata.renderers;

import xmipp.viewer.imageitems.tableitems.GalleryImageItem;
import java.awt.Component;
import java.awt.Dimension;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;

/**
 *
 * @author Juanjo Vega
 */
public class MetaDataImageRenderer extends MetaDataRowDisablerRenderer {

    public static int DEFAULT_CELL_WIDTH = 128;
    public static int DEFAULT_CELL_HEIGHT = 128;
//    private final static int BORDER_WIDTH = 5;
//    private final static int BORDER_HEIGHT = 5;
    public final static int MAX_SIZE_TO_RENDER = 512;
    private boolean renderImages = false;

    public MetaDataImageRenderer() {
        super();

        setHorizontalAlignment(JLabel.CENTER);
        setVerticalAlignment(JLabel.CENTER);
    }

    public void setRenderImages(boolean renderImages) {
        this.renderImages = renderImages;
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);

        GalleryImageItem item = (GalleryImageItem) value;
        ImageIcon icon = null;
        int w, h;

        if (renderImages && !item.isBiggerThan(MAX_SIZE_TO_RENDER)) {
            w = item.getWidth() > DEFAULT_CELL_WIDTH ? DEFAULT_CELL_WIDTH : item.getWidth();
            h = item.getWidth() > DEFAULT_CELL_HEIGHT ? DEFAULT_CELL_HEIGHT : item.getHeight();

            icon = new ImageIcon(item.getPreview(w, h).getImage());
            setPreferredSize(new Dimension(w, h));

            setText(null);

//        } else {
//            String labelStr = item.getOriginalValue();
//            String columnName = table.getColumnName(column);
//
//            int length = Math.max(labelStr.length(), columnName.length());
//
//            w = font.getSize() * length;
//            h = font.getSize() * 3 / 2;
//
//            label = getShortLabel(item.getOriginalValue(), w);
        }

        setIcon(icon);

        setToolTipText(item.getTooltipText());
        setEnabled(isRowEnabled(table, row));

        return this;
    }
//
//    private String getShortLabel(String label, int width) {
//        StringBuilder sb = new StringBuilder(label);
//        String sortLabel = sb.toString();
//
//        int w = getFontMetrics(font).stringWidth(sortLabel);
//
//        int i = 0;
//        while (w > width) {
//            sortLabel = "..." + sb.substring(i++);
//            w = getFontMetrics(font).stringWidth(sortLabel);
//        }
//
//        return sortLabel;
//    }
}
