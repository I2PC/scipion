/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package metadata.renderers;

import java.awt.Component;
import javax.swing.JLabel;
import javax.swing.JTable;
import metadata.images.TableFileItem;

/**
 *
 * @author Juanjo Vega
 */
public class MetaDataFileItemRenderer extends MetaDataRowDisablerRenderer {

    public MetaDataFileItemRenderer() {
        super();

        setHorizontalAlignment(JLabel.CENTER);
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        TableFileItem item = (TableFileItem) value;

        String str = getShortLabel(item.getOriginalValue(), table.getColumnModel().getColumn(column).getWidth());

        setToolTipText(item.getPath());

        return super.getTableCellRendererComponent(table, str, isSelected, hasFocus, row, column);
    }
}
