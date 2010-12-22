/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs.renderers;

import browser.table.micrographs.TableImageItemCOSS;
import java.awt.Component;
import javax.swing.ImageIcon;
import javax.swing.JTable;

/**
 *
 * @author Juanjo Vega
 */
public class ImageRenderer extends DisablerRenderer {

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        TableImageItemCOSS item = (TableImageItemCOSS) value;

        ImageIcon icon = new ImageIcon(item.getImagePlus().getImage());
        label.setIcon(icon);

        label.setEnabled(isRowEnabled(table, row));

        setToolTipText(item.getFileName());

        setCellBackgroundForeground(table, isSelected, hasFocus);

        return this;
    }
}
