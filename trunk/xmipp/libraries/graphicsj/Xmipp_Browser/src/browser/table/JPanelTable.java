/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * JPanelTable.java
 *
 * Created on 10-feb-2010, 13:35:14
 */
package browser.table;

import browser.ICONS_MANAGER;
import browser.imageitems.TableImageItem;
import browser.LABELS;
import browser.imageitems.FileImageItem;
import browser.imageitems.ImageItem;
import browser.imageitems.SelFileItem;
import browser.windows.ImagesWindowFactory;
import ij.ImagePlus;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.Vector;
import javax.swing.JOptionPane;
import javax.swing.ListSelectionModel;

/**
 *
 * @author Juanjo Vega
 */
public class JPanelTable extends javax.swing.JPanel {

    protected ImagesTableModel tableModel;
    protected ImagesTableRenderer tableRenderer;
    protected ImagesTableColumnModel columnModel;
    protected final static int CACHE_CONST = 2;

    /** Creates new form JPanelTable */
    public JPanelTable() {
        tableModel = new ImagesTableModel();
        columnModel = new ImagesTableColumnModel();

        initComponents();

        jtImages.setColumnModel(columnModel);
        //jtImages.setTableHeader(null);

        tableRenderer = new ImagesTableRenderer();
        jtImages.setDefaultRenderer(TableImageItem.class, tableRenderer);

        jtImages.setRowSelectionAllowed(false);
        jtImages.setColumnSelectionAllowed(false);
        jtImages.setCellSelectionEnabled(false);

        jtImages.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);

        jtImages.addMouseListener(new MouseAdapter() {

            @Override
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() > 1) {
                    openSelection();
                }
            }
        });
    }

    private void openSelection() {
        int rows[] = jtImages.getSelectedRows();
        int cols[] = jtImages.getSelectedColumns();

        for (int i = 0; i < rows.length; i++) {
            ImagesWindowFactory.openImageWindow(
                    ((TableImageItem) jtImages.getValueAt(rows[i], cols[i])).getImagePlus());
        }
    }

    protected void std_deviation() {
        Vector<TableImageItem> items = tableModel.getSelectedItems();

        if (items.size() > 0) {
            ImagePlus images[] = new ImagePlus[items.size()];

            for (int i = 0; i < items.size(); i++) {
                images[i] = items.elementAt(i).getImagePlus();
            }

            ImagePlus std_dev = ImageTableOperations.std_deviation(images);
            std_dev.setTitle("Std. Deviation");

            ImagesWindowFactory.openImageWindow(std_dev);
        } else {
            JOptionPane.showMessageDialog(this, LABELS.MESSAGE_NO_ITEMS_SELECTED);
        }
    }

    protected void mean() {
        Vector<TableImageItem> items = tableModel.getSelectedItems();

        if (items.size() > 0) {
            ImagePlus images[] = new ImagePlus[items.size()];

            // Gets the full size images.
            for (int i = 0; i < items.size(); i++) {
                images[i] = items.elementAt(i).getImagePlus();
            }

            ImagePlus mean = ImageTableOperations.mean(images);
            mean.setTitle("Mean");

            ImagesWindowFactory.openImageWindow(mean);
        } else {
            JOptionPane.showMessageDialog(this, LABELS.MESSAGE_NO_ITEMS_SELECTED);
        }
    }

    protected void normalize(boolean normalize) {
        tableModel.setNormalize(normalize);

        refreshData();
    }

    public void selectAll() {
        tableModel.selectAll();
        jtImages.repaint();
    }

    public void invertSelection() {
        tableModel.invertSelection();
        jtImages.repaint();
    }

    public void clear() {
        tableModel.clear();
        jtImages.repaint();
    }

    public void refreshData() {
        tableModel.refresh();
        jtImages.repaint();
    }

    public int getImagesCount() {
        return tableModel.getSize();
    }

    public void setZoom(int zoom) {
        tableModel.setZoom(zoom);

        refreshCacheSize();

        refreshTable();
    }

    public void setRows(int rows) {
        tableModel.setRows(rows);
    }

    public int getRows() {
        return tableModel.getRowCount();
    }

    public void setColumns(int cols) {
        tableModel.setColumns(cols);
    }

    public int getColumns() {
        return tableModel.getColumnCount();
    }

    public void addImage(ImageItem itemImage) {
        tableModel.addImage(itemImage);
    }

    public void addVolume(FileImageItem itemImage) {
        tableModel.addVolume(itemImage);
    }

    public void addVolume(SelFileItem itemImage) {
        tableModel.addVolume(itemImage);
    }

    public void setImageSelected(int index) {
        tableModel.setSelected(index);

        // Ensures item is visible.
        int indexes[] = tableModel.getRowColForIndex(index);

        jspTable.getViewport().scrollRectToVisible(jtImages.getCellRect(indexes[0], indexes[1], true));
    }

    public void setShowLabels(boolean show) {
        tableModel.setShowLabels(show);

        refreshCacheSize();

        refreshTable();
    }

    protected Dimension getCellSize() {
        int font_height = 0;
        if (tableModel.isShowingLabels()) {
            font_height = tableRenderer.getFontMetrics(tableRenderer.getFont()).getHeight();
            font_height += tableRenderer.getIconTextGap();  // Adds the extra gap.
            font_height -= jtImages.getIntercellSpacing().height;   // Removes spacing.
        }

        TableImageItem item;
        int w, h;
        int row = 0;

        do {
            item = (TableImageItem) tableModel.getItems().get(row++);

            w = item.getThumbnailWidth();
            h = item.getThumbnailHeight() + font_height;
        } while (w < 1 && h < 1 && row < tableModel.getItems().size());

        if (w < 1 || h < 1) {
            w = ICONS_MANAGER.MISSING_ITEM.getIconWidth();
            h = ICONS_MANAGER.MISSING_ITEM.getIconHeight();
            //System.out.println("[!] > w=" + w + " h=" + h);
        }

        return new Dimension(w, h);
    }

    public void refreshCacheSize() {
        if (isShowing()) {
            Dimension cellSize = getCellSize();
            Rectangle size = jspTable.getVisibleRect();

            int horizontal = (int) Math.ceil((float) size.width / (float) cellSize.width);
            int vertical = (int) Math.ceil((float) size.height / (float) cellSize.height);

            tableModel.setCacheSize(horizontal * vertical * 2);
        }
    }

    private void refreshTable() {
        if (getImagesCount() > 0) {
            Dimension dimension = getCellSize();
            tableRenderer.setPreferredSize(dimension);

            // Adjusts rows size.
            if (dimension.height > 0) {
                jtImages.setRowHeight(dimension.height);
                columnModel.setWidth(dimension.width);

                tableModel.fireTableStructureChanged();
            }

            refreshData();
        }
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jspTable = new javax.swing.JScrollPane();
        jtImages = new javax.swing.JTable();

        setLayout(new java.awt.BorderLayout());

        jtImages.setModel(tableModel);
        jtImages.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
        jtImages.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jtImagesMouseClicked(evt);
            }
        });
        jspTable.setViewportView(jtImages);

        add(jspTable, java.awt.BorderLayout.CENTER);
    }// </editor-fold>//GEN-END:initComponents

    private void jtImagesMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jtImagesMouseClicked
        int col = jtImages.columnAtPoint(evt.getPoint());
        int row = jtImages.rowAtPoint(evt.getPoint());

        // Ctrl adds items to selection, otherwise previous ones are removed.
        if (!evt.isControlDown()) {
            tableModel.clearSelection();
        }

        tableModel.toggleSelected(row, col);
        jtImages.repaint();
    }//GEN-LAST:event_jtImagesMouseClicked
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JScrollPane jspTable;
    private javax.swing.JTable jtImages;
    // End of variables declaration//GEN-END:variables
}
