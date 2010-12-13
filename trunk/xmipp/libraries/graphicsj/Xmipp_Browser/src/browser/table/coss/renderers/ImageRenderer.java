/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.coss.renderers;

import browser.table.coss.TableImageItemCOSS;
import java.awt.BorderLayout;
import java.awt.Component;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.UIManager;
import javax.swing.table.TableCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public class ImageRenderer extends JPanel implements TableCellRenderer {

    JLabel label = new JLabel();

    public ImageRenderer() {
        super();

        setLayout(new BorderLayout());
        add(label, BorderLayout.CENTER);

        label.setHorizontalAlignment(JLabel.CENTER);
        label.setVerticalAlignment(JLabel.CENTER);
    }

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        TableImageItemCOSS item = (TableImageItemCOSS) value;

        ImageIcon icon = new ImageIcon(item.getImagePlus().getImage());
        label.setIcon(icon);

//        if (path.trim().isEmpty()) {
//            label.setIcon(null);
//        } else {
//            label.setIcon(new ImageIcon(IJ.openImage(path).getImage()));//path));
//        }

        setToolTipText(item.getFileName());

        if (isSelected) {
            setBackground(table.getSelectionBackground());
            setForeground(table.getSelectionForeground());
        } else {
            setBackground(table.getBackground());
            setForeground(table.getForeground());
        }

        if (hasFocus) {
            setBorder(UIManager.getBorder("Table.focusCellHighlightBorder"));
        } else {
            setBorder(null);
        }

        return this;
    }
}
