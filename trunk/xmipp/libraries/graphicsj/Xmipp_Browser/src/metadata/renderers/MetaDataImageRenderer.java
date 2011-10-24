/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package metadata.renderers;

import browser.imageitems.tableitems.GalleryImageItem;
import ij.ImagePlus;
import java.awt.Component;
import java.awt.Font;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;

/**
 *
 * @author Juanjo Vega
 */
public class MetaDataImageRenderer extends MetaDataRowDisablerRenderer {

    private static int DEFAULT_CELL_WIDTH = 128;
    private static int DEFAULT_CELL_HEIGHT = 128;
    private final static int BORDER_WIDTH = 5;
    private final static int BORDER_HEIGHT = 5;
    private boolean renderImages;
    private Font font;

    public MetaDataImageRenderer() {
        super();

        setHorizontalAlignment(JLabel.CENTER);
        setVerticalAlignment(JLabel.CENTER);

        font = getFont();
    }

    public void setRenderImages(boolean renderImages) {
        this.renderImages = renderImages;
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        super.getTableCellRendererComponent(table, null, isSelected, hasFocus, row, column);

        GalleryImageItem item = (GalleryImageItem) value;

        String label = null;
        ImageIcon icon = null;

        if (renderImages) {
            int w = item.getWidth() > DEFAULT_CELL_WIDTH ? DEFAULT_CELL_WIDTH : item.getWidth();
            int h = item.getHeight() > DEFAULT_CELL_HEIGHT ? DEFAULT_CELL_HEIGHT : item.getHeight();

            ImagePlus preview = item.getPreview(w, h);
            icon = new ImageIcon(preview.getImage());
        } else {
            label = getShortLabel(item.getOriginalValue(), table.getColumnModel().getColumn(column).getWidth());
        }

        setText(label);
        setIcon(icon);
        setToolTipText(item.getTooltipText());
        setEnabled(isRowEnabled(table, row));

        return this;
    }

    public int getCellWidth() {
        return DEFAULT_CELL_WIDTH;
    }

    public int getCellHeight() {
        int height;

        if (renderImages) {
            height = DEFAULT_CELL_HEIGHT;
        } else {
            height = font.getSize();
            height += height * 0.5; // Extra gap.
        }

        return height;
    }

    private String getShortLabel(String label, int width) {
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
}
