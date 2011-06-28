/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.renderers;

import browser.imageitems.tableitems.AbstractTableImageItem;
import ij.ImagePlus;
import java.awt.Component;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;
import browser.table.models.AbstractXmippTableModel;
import java.awt.Color;
import java.awt.Font;
import java.awt.Image;
import javax.swing.BorderFactory;
import javax.swing.border.Border;
import javax.swing.table.DefaultTableCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public class ImageRenderer extends DefaultTableCellRenderer {

    protected final static String COLOR_FILL = "#3aff35";
    protected final static String COLOR_OUTLINE = "#2aa430";
    protected Border BORDER_SELECTED = BorderFactory.createLineBorder(Color.RED, 1);
    protected Border BORDER_FOCUSED = BorderFactory.createLineBorder(Color.RED, 3);

//    protected Image image; <- Check if will use it again: Defined below.
//    protected Polygon triangleMark;
    @Override
    public Component getTableCellRendererComponent(JTable table, Object object, boolean isSelected, boolean hasFocus, int row, int column) {
        AbstractTableImageItem item = (AbstractTableImageItem) object;

        // Calls super class so foreground, background, borders and rest of stuff is set.
        super.getTableCellRendererComponent(table, null,
                item != null && item.isSelected(),
                item != null && hasFocus, row, column);

        if (item != null) {
//        DEBUG.printMessage("*** Rendering: " + item + " S: " + "<item.slice>" + " N: " + item.nimage);
            AbstractXmippTableModel tableModel = (AbstractXmippTableModel) table.getModel();

            // Loads image...
            ImagePlus imp = item.getPreview();

            // ... and sets it.
            setEnabled(item.isEnabled());

            // Normalizes image (if sets in tablemodel)
            normalize(imp, tableModel);

            Image image = imp.getImage();
            setIcon(new ImageIcon(image));

            setOpaque(true);
            setHorizontalAlignment(JLabel.CENTER);
            setHorizontalTextPosition(JLabel.CENTER);
            setVerticalTextPosition(JLabel.BOTTOM);

            // Tooltip.
            setToolTipText(item.getTooltipText());

            // (Shows label only when required).
            if (tableModel.isShowingLabels()) {
                String label = cutString(
                        String.valueOf(item.getLabelValue(tableModel.getSelectedLabel())),
                        table.getColumnModel().getColumn(column).getWidth());

                setText(label);
            } else {
                setText(null);
            }

            // Hacking borders to enhance the default one.
            if (item.isSelected()) {
                setBorder(BORDER_SELECTED);
            }

            if (hasFocus) {
                setBorder(BORDER_FOCUSED);
            }
        } else {
            setIcon(null);
            setText(null);
            setToolTipText(null);
        }

        return this;
    }

    private void normalize(ImagePlus image, AbstractXmippTableModel tableModel) {
        if (tableModel.isNormalizing()) {
            image.getProcessor().setMinAndMax(tableModel.getNormalizeMin(),
                    tableModel.getNormalizeMax());
        } else {
            image.getProcessor().resetMinAndMax();
        }

        image.updateImage();  // Repaint
    }

    protected String cutString(String string, int columnWidth) {
        StringBuffer sb = new StringBuffer(string);
        String str = sb.toString();

        Font font = getFont();
        int w = getFontMetrics(font).stringWidth(str);

        int i = 0;
        while (w > columnWidth) {
            str = "..." + sb.substring(i++);
            w = getFontMetrics(font).stringWidth(str);
        }

        return str;
    }
    /*
    @Override
    public void paint(Graphics g) {
    super.paint(g);

    boolean mark = Math.random() < 0.5;

    if (mark) {
    //if (triangleMark == null) {
    triangleMark = createPolygon(getWidth(), getHeight());
    //}

    // Area and...
    g.setColor(Color.decode(COLOR_FILL));
    g.fillPolygon(triangleMark);

    // ...outline.
    g.setColor(Color.decode(COLOR_OUTLINE));
    g.drawPolygon(triangleMark);
    }
    }

    public static Polygon createPolygon(int w, int h) {
    Polygon t = new Polygon();

    t.addPoint(w - w / 5, h - 1);
    t.addPoint(w - 1, h - h / 5);
    t.addPoint(w - 1, h - 1);

    return t;
    }*/
}
