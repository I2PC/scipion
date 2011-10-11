/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package metadata.renderers;

import java.awt.Component;
import java.text.DecimalFormat;
import javax.swing.JLabel;
import javax.swing.JTable;

/**
 *
 * @author Juanjo Vega
 */
public class MetaDataDoubleRenderer extends MetaDataRowDisablerRenderer {

    protected final static DecimalFormat formatter = new DecimalFormat("#.##E0");

    public MetaDataDoubleRenderer() {
        super();

        setHorizontalAlignment(JLabel.RIGHT);
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        return super.getTableCellRendererComponent(table, "<b>" + formatter.format(value).replace('E', 'e') + "</b>", isSelected, hasFocus, row, column);
    }
}
