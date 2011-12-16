/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.gallery.renderers;

import browser.imageitems.tableitems.AbstractGalleryImageItem;
import ij.ImagePlus;
import java.awt.Component;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;
import browser.gallery.models.AbstractXmippTableModel;
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

    Border BORDER_SELECTED = new StrokeBorder(Color.RED, 3);
    Border BORDER_FOCUSED = BorderFactory.createLineBorder(Color.RED, 3);

    public ImageRenderer() {
        super();
        setOpaque(true);
        setHorizontalAlignment(JLabel.CENTER);
        setHorizontalTextPosition(JLabel.CENTER);
        setVerticalTextPosition(JLabel.BOTTOM);
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object object, boolean isSelected, boolean hasFocus, int row, int column) {
        AbstractGalleryImageItem item = (AbstractGalleryImageItem) object;

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
            try{
            normalize(imp, tableModel);
            }catch(Exception ex){
                System.out.println(ex.getMessage());
            }
            Image image = imp.getImage();
            setIcon(new ImageIcon(image));

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
        StringBuilder sb = new StringBuilder(string);
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
}
