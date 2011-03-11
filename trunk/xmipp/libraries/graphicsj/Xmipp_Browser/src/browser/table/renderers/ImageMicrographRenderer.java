/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.renderers;

import browser.imageitems.TableImageItem;
import ij.ImagePlus;
import java.awt.Component;
import javax.swing.ImageIcon;
import javax.swing.JTable;

/**
 *
 * @author Juanjo Vega
 */
public class ImageMicrographRenderer extends RowDisablerRenderer {

    public static final int CELL_WIDTH = 128;
    public static final int CELL_HEIGHT = 128;

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        TableImageItem item = (TableImageItem) value;

        ImagePlus preview = item.getPreview(CELL_WIDTH, CELL_HEIGHT);

        ImageIcon icon = new ImageIcon(preview.getImage());
        setIcon(icon);

        setEnabled(isRowEnabled(table, row));

        setToolTipText(item.getFileName());

        setCellBackgroundForeground(table, isSelected, hasFocus);

        return this;
    }
}
