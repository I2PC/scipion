/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package rotspectra;

import java.awt.Color;
import java.awt.Component;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.UIManager;
import javax.swing.table.DefaultTableCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public class RotSpectraRenderer extends DefaultTableCellRenderer {

    int imagew, imageh;

    public RotSpectraRenderer(int imagew, int imageh) {
        this.imagew = imagew;
        this.imageh = imageh;
        setSize(imagew, imageh);

        setBackground(Color.WHITE);
        setOpaque(true);
        //setHorizontalAlignment(JLabel.CENTER);
        setHorizontalTextPosition(JLabel.CENTER);
        setVerticalTextPosition(JLabel.BOTTOM);
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object o, boolean isSelected, boolean hasFocus, int rowIndex, int colIndex) {
//        super.getTableCellRendererComponent(table, o, bln, bln1, i, i1);

//        setBorder(isSelected ? UIManager.getBorder("List.focusCellHighlightBorder") : null);
        setBorder(hasFocus ? UIManager.getBorder("List.focusCellHighlightBorder") : null);

        // Sets plot image.
        RotSpectraVector item = (RotSpectraVector) o;
        setIcon(new ImageIcon(item.getPreview(imagew, imageh)));

        RotSpectraTableModel tableModel = (RotSpectraTableModel) table.getModel();
        if (tableModel.isShowingLabels()) {
            setText("<html><b>"+item.getNImages() + "</b> images</html>");
        } else {
            setText(null);
        }

        // Tooltip.
        setToolTipText(item.getTooltipText());

        return this;
    }
}
