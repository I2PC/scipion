/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package rotspectra;

import java.awt.Component;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public class RotSpectraRenderer extends DefaultTableCellRenderer {

    int width, height;

    public RotSpectraRenderer(int width, int height) {
        this.width = width;
        this.height = height;
        setSize(width, height);
    }

    @Override
    public int getWidth() {
        return width;
    }

    @Override
    public int getHeight() {
        return height;
    }

    @Override
    public Component getTableCellRendererComponent(JTable jtable, Object o, boolean bln, boolean bln1, int i, int i1) {
        super.getTableCellRendererComponent(jtable, o, bln, bln1, i, i1);

        // Sets plot image.
        RotSpectraVector item = (RotSpectraVector) o;
        setIcon(new ImageIcon(item.getPlotImage(width, height)));
        setText(null);

        return this;
    }
}
