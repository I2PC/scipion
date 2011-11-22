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
import metadata.models.MetaDataTableModel;
import xmipp.MDLabel;

/**
 *
 * @author Juanjo Vega
 */
public class MetaDataImageRenderer extends MetaDataRowDisablerRenderer {

//    private static int DEFAULT_CELL_WIDTH = 128;
//    private static int DEFAULT_CELL_HEIGHT = 128;
//    private final static int BORDER_WIDTH = 5;
//    private final static int BORDER_HEIGHT = 5;
    private boolean renderImages;
    MetaDataTableModel tableModel;
    private Font font;

    public MetaDataImageRenderer(MetaDataTableModel tableModel) {
        super();

        this.tableModel = tableModel;

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

        if (renderImages && item.getLabel() != MDLabel.MDL_IMAGE) {   // Image is not rendered.
            ImagePlus preview = item.getPreview(getCellWidth(), getCellHeight());
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

    public int getCellWidth(String columnTitle) {
        int width = getCellWidth();

        if (width < 0) {
            width = font.getSize() * columnTitle.length();
            width += width * 0.5; // Extra gap.
        }

        return width;
    }

    int getCellWidth() {
        return renderImages ? tableModel.getMaxColumnWidth() : -1;
    }

    public int getCellHeight() {
        int height;

        if (renderImages) {
            height = tableModel.getMaxRowHeight();
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
