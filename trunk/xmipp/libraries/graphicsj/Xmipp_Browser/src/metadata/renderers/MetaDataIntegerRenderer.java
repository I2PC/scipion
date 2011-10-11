/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package metadata.renderers;

import java.awt.Component;
import javax.swing.JLabel;
import javax.swing.JTable;

/**
 *
 * @author Juanjo Vega
 */
public class MetaDataIntegerRenderer extends MetaDataRowDisablerRenderer {

//    protected final static DecimalFormat formatter = new DecimalFormat("#.##E0");
    public MetaDataIntegerRenderer() {
        super();

        setHorizontalAlignment(JLabel.RIGHT);
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        return super.getTableCellRendererComponent(table, "<b>" + (Integer) value + "</b>", isSelected, hasFocus, row, column);
    }
}
