/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */

/*
 * JPanelTable.java
 *
 * Created on 10-feb-2010, 13:35:14
 */
package browser.table;

import browser.imageitems.TableImageItem;
import browser.LABELS;
import browser.imageitems.listitems.AbstractImageItem;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.util.Vector;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.JSeparator;
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;

/**
 *
 * @author Juanjo Vega
 */
public class JPanelTable extends javax.swing.JPanel {

    protected final static int PIXEL_SIZE = 4;    // 4B/px
    protected final static int CACHE_MAX_SIZE = 200 * 1024 * 1024;    // 1MB = 1024 KB => 200MB * 1024^2
    protected ImagesTableModel tableModel;
    protected ImagesTableRenderer tableRenderer;
    protected ImagesTableColumnModel columnModel;
    protected JPopupMenu jPopupMenu = new JPopupMenu();
    protected JMenuItem jmiEnable = new JMenuItem(LABELS.LABEL_TABLE_ENABLE);
    protected JMenuItem jmiDisable = new JMenuItem(LABELS.LABEL_TABLE_DISABLE);
    protected JMenuItem jmiEnableAll = new JMenuItem(LABELS.LABEL_TABLE_ENABLE_ALL);
    protected JMenuItem jmiDisableAll = new JMenuItem(LABELS.LABEL_TABLE_DISABLE_ALL);
    protected JMenuItem jmiSaveAsImages = new JMenuItem(LABELS.LABEL_TABLE_SAVE_AS_IMAGES);
    protected JMenuItem jmiSaveAsStack = new JMenuItem(LABELS.LABEL_TABLE_SAVE_AS_STACK);
    protected JMenuItem jmiSaveAsSelfile = new JMenuItem(LABELS.LABEL_TABLE_SAVE_AS_SELFILE);
    private Vector<TableImageItem> selectedItems;    // Currently selected items (used when right-click)
    private boolean autoAdjustColumns = false;

    /** Creates new form JPanelTable */
    public JPanelTable() {
        tableModel = new ImagesTableModel();
        columnModel = new ImagesTableColumnModel();

        initComponents();

        jtImages.setColumnModel(columnModel);
        jtImages.setIntercellSpacing(new Dimension(2, 2));

        tableRenderer = new ImagesTableRenderer();
        jtImages.setDefaultRenderer(TableImageItem.class, tableRenderer);

        jtImages.setRowSelectionAllowed(false);
        jtImages.setColumnSelectionAllowed(false);
        jtImages.setCellSelectionEnabled(false);

        jtImages.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);

        jPopupMenu.add(jmiEnable);
        jPopupMenu.add(jmiDisable);
        jPopupMenu.add(new JSeparator());
        jPopupMenu.add(jmiDisableAll);
        jPopupMenu.add(jmiEnableAll);
        jPopupMenu.add(new JSeparator());
        jPopupMenu.add(new JSeparator());
        jPopupMenu.add(jmiSaveAsImages);
        jPopupMenu.add(jmiSaveAsStack);
        jPopupMenu.add(jmiSaveAsSelfile);

        jmiEnable.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, InputEvent.CTRL_DOWN_MASK));
        jmiDisable.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_D, InputEvent.CTRL_DOWN_MASK));

        jtImages.addKeyListener(new KeyListener() {

            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                if (e.isControlDown()) {
                    if (e.getKeyCode() == KeyEvent.VK_D) {
                        disableSelectedItems();
                    }
                    if (e.getKeyCode() == KeyEvent.VK_E) {
                        enableSelectedItems();
                    }
                }
            }
        });

        jmiEnable.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                enableSelectedItems();
            }
        });

        jmiDisable.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                disableSelectedItems();
            }
        });

        jmiEnableAll.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                enableAllItems();
            }
        });

        jmiDisableAll.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                disableAllItems();
            }
        });

        jmiSaveAsImages.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                if (!saveAsImages()) {
                    JOptionPane.showMessageDialog(getParent(), LABELS.MESSAGE_NO_ITEMS_SELECTED);
                }
            }
        });

        jmiSaveAsStack.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                if (!saveAsStack()) {
                    JOptionPane.showMessageDialog(getParent(), LABELS.MESSAGE_NO_ITEMS_SELECTED);
                }
            }
        });

        jmiSaveAsSelfile.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                if (!saveAsSelFile()) {
                    JOptionPane.showMessageDialog(getParent(), LABELS.MESSAGE_NO_ITEMS_SELECTED);
                }
            }
        });
    }

    protected int getDisplayableColumns() {
        return getVisibleRect().width / getCellWidth();
    }

    private boolean saveAsSelFile() {
        selectedItems = tableModel.getSelectedAndEnabledItems();

        if (selectedItems != null) {
            if (selectedItems.size() > 0) {
                int w = selectedItems.elementAt(0).getImagePlus().getWidth();
                int h = selectedItems.elementAt(0).getImagePlus().getHeight();
                ImageStack stack = new ImageStack(w, h, selectedItems.size());

                for (int i = 0; i < selectedItems.size(); i++) {
                    ImagePlus ip = selectedItems.elementAt(i).getImagePlus();
                    stack.setPixels(ip.getProcessor().getPixels(), i + 1);
                }

                ImagePlus ip = new ImagePlus("", stack);
                IJ.run(ip, "Sel writer", "");

                return true;
            }
        }

        return false;
    }

    private boolean saveAsStack() {
        selectedItems = tableModel.getSelectedAndEnabledItems();

        if (selectedItems != null) {
            if (selectedItems.size() > 0) {
                int w = selectedItems.elementAt(0).getImagePlus().getWidth();
                int h = selectedItems.elementAt(0).getImagePlus().getHeight();
                ImageStack stack = new ImageStack(w, h, selectedItems.size());

                for (int i = 0; i < selectedItems.size(); i++) {
                    ImagePlus ip = selectedItems.elementAt(i).getImagePlus();
                    stack.setPixels(ip.getProcessor().getPixels(), i + 1);
                }

                ImagePlus ip = new ImagePlus("", stack);
                IJ.run(ip, "Spider writer", "");

                return true;
            }
        }

        return false;
    }

    private boolean saveAsImages() {
        selectedItems = tableModel.getSelectedAndEnabledItems();

        if (selectedItems.size() > 0) {
            for (int i = 0; i < selectedItems.size(); i++) {
                IJ.run(selectedItems.elementAt(i).getImagePlus(), "Spider writer", "");
            }

            return true;
        }

        return false;
    }

    private void enableItems(boolean enable, Vector<TableImageItem> items) {
        if (items != null) {
            for (int i = 0; i < items.size(); i++) {
                items.elementAt(i).setEnabled(enable);
            }

            jtImages.repaint();
        }
    }

    private void enableAllItems() {
        enableItems(true, tableModel.getItems());
    }

    private void disableAllItems() {
        enableItems(false, tableModel.getItems());
    }

    private void enableSelectedItems() {
        enableItems(true, selectedItems);
    }

    private void disableSelectedItems() {
        enableItems(false, selectedItems);
    }

    private void openSelection() {
        int rows[] = jtImages.getSelectedRows();
        int cols[] = jtImages.getSelectedColumns();

        for (int i = 0; i < rows.length; i++) {
            ImagesWindowFactory.openImageWindow(
                    (TableImageItem) jtImages.getValueAt(rows[i], cols[i]));
        }
    }

    protected void average() {
        Vector<TableImageItem> items = tableModel.getSelectedAndEnabledItems();

        if (items.size() > 0) {
            ImagePlus average = ImageTableOperations.average(items);
            average.setTitle("Average");

            ImagesWindowFactory.openImageWindow(average, false);
        } else {
            JOptionPane.showMessageDialog(this, LABELS.MESSAGE_NO_ITEMS_SELECTED);
        }
    }

    protected void std_deviation() {
        Vector<TableImageItem> items = tableModel.getSelectedAndEnabledItems();

        if (items.size() > 0) {
            ImagePlus std_dev = ImageTableOperations.std_deviation(items);
            std_dev.setTitle("Std. Deviation");

            ImagesWindowFactory.openImageWindow(std_dev, false);
        } else {
            JOptionPane.showMessageDialog(this, LABELS.MESSAGE_NO_ITEMS_SELECTED);
        }
    }

    public void normalize(double min, double max) {
        tableModel.setNormalized(min, max);

        refreshData();
    }

    protected void setNormalizedAuto() {
        tableModel.setNormalizedAuto();

        refreshData();
    }

    public void disableNormalization() {
        tableModel.disableNormalization();

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

    public int getItemsCount() {
        return tableModel.getSize();
    }

    public void setZoom(int zoom) {
        tableModel.setZoom(zoom);

        refreshCacheSize();

        refreshTable();
    }

    public void setAutoAdjustColumns(boolean autoAdjustColumns) {
        this.autoAdjustColumns = autoAdjustColumns;

        autoAdjustColumns();
    }

    public Vector<TableImageItem> getItems() {
        return tableModel.getItems();
    }

    public void setTableSize(int rows, int cols) {
        tableModel.setTableSize(rows, cols);
    }

    public void setRows(int rows) {
        tableModel.setRows(rows);
    }

    public int getRowsCount() {
        return tableModel.getRowCount();
    }

    public void setColumns(int cols) {
        tableModel.setColumns(cols);
    }

    public int getColumnsCount() {
        return tableModel.getColumnCount();
    }

    public void addImage(AbstractImageItem itemImage) {
        tableModel.addImage(itemImage);
    }

    public void addVolume(AbstractImageItem itemImage) {
        tableModel.addVolume(itemImage);
    }

    public int addVolumeFile(File file) {
        return tableModel.addVolumeFile(file);
    }

    public void addSelVolume(String dir, String file) {
        tableModel.addSelVolume(dir, file);
    }

    public int getCellWidth() {
        return tableModel.getItems().elementAt(0).getThumbnailWidth();//ImagePlus().getHeight();
    }

    public int getCellHeight() {
        return tableModel.getItems().elementAt(0).getThumbnailHeight();//ImagePlus().getHeight();
    }

    public void setImageSelected(int index) {
        tableModel.setSelected(index);

        int indexes[] = tableModel.getRowColForIndex(index);

        // Gets current selected cell bounds.
        Rectangle rect = jtImages.getCellRect(indexes[0], indexes[1], true);

        // Ensures item is visible.
        Point pos = jspTable.getViewport().getViewPosition();
        rect.translate(-pos.x, -pos.y);
        jspTable.getViewport().scrollRectToVisible(rect);

        repaint();
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
        } while (row < tableModel.getItems().size());

        return new Dimension(w, h);
    }

    public void refreshCacheSize() {
        if (isShowing()) {
            Dimension cellSize = getCellSize();
//            Rectangle size = jspTable.getVisibleRect();

//            int horizontal = (int) Math.ceil((float) size.width / (float) cellSize.width);
//            int vertical = (int) Math.ceil((float) size.height / (float) cellSize.height);

            // new_cache_size = (200 MB = 200 * 1024 * 8) / (cell.w * cell.h * 8)

            double siCell = cellSize.width * cellSize.height * PIXEL_SIZE;
            double cacheSize = CACHE_MAX_SIZE / siCell;
//            System.out.println(" >> Cell W: " + cellSize.width + " / H:" + cellSize.height);
//            System.out.println(" >> Cell size: " + siCell / 1024 + " KB");
//            System.out.println(" >> Cache size: " + CACHE_MAX_SIZE / (1024 * 1024) + " MB -> " + CACHE_MAX_SIZE / 1024 + " KB");
//            System.out.println(" >> Max items in cache: " + cacheSize);
//            System.out.println(" ##############################################################");

            //tableModel.setCacheSize(horizontal * vertical * 2);
            tableModel.resizeCache((int) cacheSize);
        }
    }

    private void refreshTable() {
        if (getItemsCount() > 0) {
            Dimension dimension = getCellSize();
            tableRenderer.setPreferredSize(dimension);

            // Adjusts rows size.
            jtImages.setRowHeight(dimension.height > 0 ? dimension.height : 1);
            columnModel.setWidth(dimension.width > 0 ? dimension.width : 1);

            tableModel.fireTableStructureChanged();

            refreshData();

            autoAdjustColumns();
        }
    }

    protected void autoAdjustColumns() {
        if (autoAdjustColumns) {
            int displayableColumns = getDisplayableColumns();

            if (getColumnsCount() != displayableColumns) {
                setColumns(displayableColumns);
            }
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
        if (evt.getButton() == MouseEvent.BUTTON1) {  // Left click.
            if (evt.getClickCount() > 1) {
                openSelection();
            } else {
                int col = jtImages.columnAtPoint(evt.getPoint());
                int row = jtImages.rowAtPoint(evt.getPoint());

                // Ctrl adds items to selection, otherwise previous ones are removed.
                if (!evt.isControlDown()) {
                    tableModel.clearSelection();
                }

                // Right click selects an item (but doesn't deselect it)
//            if (evt.getButton() == MouseEvent.BUTTON3 && !((TableImageItem) tableModel.getValueAt(row, col)).selected) {
                tableModel.toggleSelected(row, col);
                //tableModel.setSelected(row, col, true);
//            }
//            ((TableImageItem) tableModel.getValueAt(row, col)).selected = evt.getButton() == MouseEvent.BUTTON3;

                jtImages.repaint();
            }
        } else if (evt.getButton() == MouseEvent.BUTTON3) {  // Right click.
            selectedItems = tableModel.getSelectedItems();

            if (selectedItems != null) {
                jPopupMenu.show(evt.getComponent(), evt.getX(), evt.getY());
            }
        }

        /*
        if (evt.getButton() == MouseEvent.BUTTON1) {
        int col = jtImages.columnAtPoint(evt.getPoint());
        int row = jtImages.rowAtPoint(evt.getPoint());

        // Ctrl adds items to selection, otherwise previous ones are removed.
        if (!evt.isControlDown()) {
        tableModel.clearSelection();
        }

        // Right click selects an item (but doesn't deselect it)
        //            if (evt.getButton() == MouseEvent.BUTTON3 && !((TableImageItem) tableModel.getValueAt(row, col)).selected) {
        tableModel.toggleSelected(row, col);
        //tableModel.setSelected(row, col, true);
        //            }
        //            ((TableImageItem) tableModel.getValueAt(row, col)).selected = evt.getButton() == MouseEvent.BUTTON3;

        jtImages.repaint();
        }*/
    }//GEN-LAST:event_jtImagesMouseClicked
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JScrollPane jspTable;
    private javax.swing.JTable jtImages;
    // End of variables declaration//GEN-END:variables
}
