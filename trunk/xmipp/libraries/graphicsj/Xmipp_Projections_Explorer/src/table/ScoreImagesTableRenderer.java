/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package table;

import java.awt.Color;
import java.awt.Component;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.border.Border;
import javax.swing.table.TableCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public class ScoreImagesTableRenderer extends JLabel implements TableCellRenderer {

    //private Border borderSelected = BorderFactory.createLineBorder(Color.RED, 20);
    //private Border borderFocused = BorderFactory.createLineBorder(Color.RED, 3);
    private Border bordergood = BorderFactory.createLineBorder(Color.BLUE, 1);
    private Border borderbad = BorderFactory.createLineBorder(Color.RED, 1);

    public ScoreImagesTableRenderer() {
        super();
    }

    public Component getTableCellRendererComponent(JTable table, Object object, boolean isSelected, boolean hasFocus, int row, int column) {
        ScoreItem item = (ScoreItem) object;

        if (item != null) {
            setIcon(new ImageIcon(item.getImage()));

            setOpaque(true);
            setHorizontalAlignment(CENTER);
            setHorizontalTextPosition(CENTER);
            setVerticalTextPosition(BOTTOM);

            // Tooltip.
            setToolTipText(item.getTooltip());
        }

        // (Shows label only when required).
        setText(item.getLabel());
        /*
        if (isSelected) {    // If is selected...
        System.out.println(" *** selected: " + item.getLabel());
        if (hasFocus) {    // ...and focused as well.
        setBorder(borderFocused);
        setBackground(table.getSelectionBackground());
        setForeground(table.getSelectionForeground());
        } else {    // ...otherwise.
        setBorder(borderSelected);
        setBackground(table.getSelectionBackground());
        setForeground(table.getSelectionForeground());
        }
        } else {*/
        setBorder(item.good ? bordergood : borderbad);
        setBackground(table.getBackground());
        setForeground(table.getForeground());
//        }

        return this;
    }
}
