/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.renderers;

import browser.ICONS_MANAGER;
import browser.imageitems.TableImageItem;
import browser.table.ImagesTableModel;
import ij.ImagePlus;
import java.awt.Color;
import java.awt.Component;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.UIManager;
import javax.swing.border.BevelBorder;
import javax.swing.border.Border;

/**
 *
 * @author Juanjo Vega
 */
public class ImageRenderer extends ImageMicrographRenderer {

    protected final static Border BORDER_ENABLED = BorderFactory.createBevelBorder(BevelBorder.RAISED);
    protected final static Border BORDER_DISABLED = null;//BorderFactory.createBevelBorder(EtchedBorder.LOWERED);
    protected final static Border BORDER_SELECTED = BorderFactory.createLineBorder(Color.RED, 1);
    protected final static Border BORDER_FOCUSED = BorderFactory.createLineBorder(Color.RED, 3);
    protected final static Border BORDER_ENABLED_SELECTED = BorderFactory.createCompoundBorder(BORDER_ENABLED, BORDER_SELECTED);
    protected final static Border BORDER_ENABLED_FOCUSED = BorderFactory.createCompoundBorder(BORDER_ENABLED, BORDER_FOCUSED);
    protected final static Border BORDER_DISABLED_SELECTED = BorderFactory.createCompoundBorder(BORDER_DISABLED, BORDER_SELECTED);
    protected final static Border BORDER_DISABLED_FOCUSED = BorderFactory.createCompoundBorder(BORDER_DISABLED, BORDER_FOCUSED);
    protected final static Color BACKGROUND = UIManager.getColor("Table.background");
    protected final static Color BACKGROUND_SELECTED = UIManager.getColor("Table.selectionBackground");
    protected final static Color FOREGROUND = UIManager.getColor("Table.foreground");
    protected final static Color FOREGROUND_SELECTED = UIManager.getColor("Table.selectionForeground");

    @Override
    public Component getTableCellRendererComponent(JTable table, Object object, boolean isSelected, boolean hasFocus, int row, int column) {
        TableImageItem item = (TableImageItem) object;

        if (item != null) {
            //System.out.println("*** Rendering: " + item + " S: " + item.slice + " N: " + item.nimage);

            // Loads image...
            ImagePlus img = item.getPreview();

            // ... and sets it.
            if (img != null) {
                setEnabled(item.isEnabled());

                setIcon(new ImageIcon(img.getImage()));
            } else {
                setIcon(ICONS_MANAGER.MISSING_ITEM);
            }

            setOpaque(true);
            setHorizontalAlignment(JLabel.CENTER);
            setHorizontalTextPosition(JLabel.CENTER);
            setVerticalTextPosition(JLabel.BOTTOM);

            // Tooltip.
            setToolTipText(item.getTooltipText());

            // (Shows label only when required).
            if (((ImagesTableModel) table.getModel()).isShowingLabels()) {
                setText(item.getLabel());
            } else {
                setText(null);
            }

            if (item.isEnabled()) {
                if (hasFocus) {
                    setBorder(BORDER_ENABLED_FOCUSED);
                } else if (item.isSelected()) {
                    setBorder(BORDER_ENABLED_SELECTED);
                } else {
                    setBorder(BORDER_ENABLED);
                }
            } else {
                if (hasFocus) {
                    setBorder(BORDER_DISABLED_FOCUSED);
                } else if (item.isSelected()) {
                    setBorder(BORDER_DISABLED_SELECTED);
                } else {
                    setBorder(BORDER_DISABLED);
                }
            }
        } else {
            setIcon(null);
            setText(null);
            setBorder(null);
            setBackground(BACKGROUND);
            setForeground(FOREGROUND);
        }

        return this;
    }
}
