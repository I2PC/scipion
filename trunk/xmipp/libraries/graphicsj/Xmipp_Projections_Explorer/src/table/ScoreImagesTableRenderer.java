/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package table;

import java.awt.Color;
import java.awt.Component;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JTable;
import javax.swing.border.Border;
import javax.swing.table.DefaultTableCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public class ScoreImagesTableRenderer extends DefaultTableCellRenderer {//JLabel implements TableCellRenderer {

    private Border bordergood = BorderFactory.createLineBorder(Color.BLUE, 1);
    private Border borderbad = BorderFactory.createLineBorder(Color.RED, 1);

    @Override
    public Component getTableCellRendererComponent(JTable table, Object object, boolean isSelected, boolean hasFocus, int row, int column) {
        super.getTableCellRendererComponent(table, object, isSelected, hasFocus, row, column);
        
        ScoreItem item = (ScoreItem) object;

        if (item != null) {
            setIcon(new ImageIcon(item.getImagePlus().getImage()));
            

            setOpaque(true);
            setHorizontalAlignment(CENTER);
            setHorizontalTextPosition(CENTER);
            setVerticalTextPosition(BOTTOM);

            // Tooltip.
            setToolTipText(item.getTooltip());
        }

        // (Shows label only when required).
        setText(item.getLabel());

        setBorder(item.good ? bordergood : borderbad);
        setBackground(table.getBackground());
        setForeground(table.getForeground());

        return this;
    }
}
