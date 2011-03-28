/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package table;

import constants.LABELS;
import java.awt.Color;
import java.awt.Component;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;
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
