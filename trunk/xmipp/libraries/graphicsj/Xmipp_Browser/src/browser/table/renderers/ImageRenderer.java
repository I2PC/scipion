/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.renderers;

import browser.ICONS_MANAGER;
import browser.imageitems.TableImageItem;
import ij.ImagePlus;
import java.awt.Component;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;
import browser.table.ImagesTableModel;
import java.awt.Color;
import javax.swing.BorderFactory;
import javax.swing.border.Border;
import javax.swing.table.DefaultTableCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public class ImageRenderer extends DefaultTableCellRenderer {

    protected Border BORDER_SELECTED = BorderFactory.createLineBorder(Color.RED, 1);
    protected Border BORDER_FOCUSED = BorderFactory.createLineBorder(Color.RED, 3);
    protected boolean showLabels = false;

    public void setShowLabels(boolean showLabels) {
        this.showLabels = showLabels;
    }

    public boolean isShowingLabels() {
        return showLabels;
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object object, boolean isSelected, boolean hasFocus, int row, int column) {
        TableImageItem item = (TableImageItem) object;

        // Calls super class so foreground, background, borders and rest of stuff is set.
        super.getTableCellRendererComponent(table, null,
                item != null && item.isSelected(),
                item != null && hasFocus, row, column);

        if (item != null) {
            //System.out.println("*** Rendering: " + item + " S: " + item.slice + " N: " + item.nimage);
            ImagesTableModel tableModel = (ImagesTableModel) table.getModel();

            // Loads image...
            ImagePlus img = item.getPreview();//tableModel.getZoomScale());

            // ... and sets it.
            if (img != null) {
                setEnabled(item.isEnabled());

                // Normalizes image (if sets in tablemodel)
                normalize(img, tableModel);

                setIcon(new ImageIcon(img.getImage()));
            } else {
                setIcon(ICONS_MANAGER.MISSING_ITEM);
            }

            setOpaque(true);
            setHorizontalAlignment(JLabel.CENTER);
            setHorizontalTextPosition(JLabel.CENTER);
            setVerticalTextPosition(JLabel.BOTTOM);

            // Tooltip.
            setToolTipText(item.getLabel());//TooltipText());

            // (Shows label only when required).
            if (isShowingLabels()) {
                setText(item.getLabel());
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
        }

        return this;
    }

    private void normalize(ImagePlus image, ImagesTableModel tableModel) {
        if (tableModel.isNormalizing()) {
            image.getProcessor().setMinAndMax(tableModel.getNormalizeMin(),
                    tableModel.getNormalizeMax());
        } else {
            image.getProcessor().resetMinAndMax();
        }

        image.updateImage();  // Repaint
    }
}
