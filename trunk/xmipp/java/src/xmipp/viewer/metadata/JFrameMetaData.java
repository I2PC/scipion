/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * JFrameMetaData.java
 *
 * Created on Oct 25, 2010, 10:50:01 AM
 */
package xmipp.viewer.metadata;

import ij.IJ;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.ItemEvent;
import java.io.File;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JTable;
import javax.swing.LookAndFeel;
import javax.swing.SwingUtilities;
import javax.swing.event.ListSelectionEvent;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableRowSorter;
import xmipp.utils.DEBUG;
import xmipp.utils.Labels;
import xmipp.viewer.imageitems.tableitems.GalleryImageItem;
import xmipp.viewer.windows.ImagesWindowFactory;
import java.awt.Component;
import java.awt.event.KeyListener;
import javax.swing.JTextField;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import xmipp.viewer.metadata.renderers.MetaDataStringRenderer;
import xmipp.viewer.metadata.renderers.MetaDataIntegerRenderer;
import xmipp.viewer.metadata.renderers.MetaDataImageRenderer;
import xmipp.viewer.metadata.renderers.MetaDataFileItemRenderer;
import xmipp.viewer.metadata.renderers.MetaDataDoubleRenderer;
import xmipp.viewer.metadata.renderers.RowHeaderRenderer;
import xmipp.viewer.metadata.renderers.editors.TableFileItemEditor;
import xmipp.viewer.metadata.filters.RowEnableFilter;
import xmipp.viewer.metadata.models.MetaDataTableModel;
import xmipp.viewer.metadata.models.XTableColumnModel;
import xmipp.viewer.metadata.models.ImagesRowHeaderModel;
import xmipp.viewer.metadata.dialogs.JDialogColumnsSelector;
import xmipp.viewer.metadata.dialogs.JDialogTextFile;
import xmipp.viewer.metadata.images.TableFileItem;
import xmipp.viewer.metadata.images.TableMetaDataItem;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameMetaData extends JFrame {

    private JTable table;
    private MetaDataTableModel tableModel;
    private ImagesRowHeaderModel rowHeaderModel;
    private XTableColumnModel columnModel = new XTableColumnModel();
    private boolean canRender;
    private JFileChooser fc;
    private JList rowHeader;
    private RowHeaderRenderer rowHeaderRenderer = new RowHeaderRenderer();
    private MetaDataFileItemRenderer fileRenderer = new MetaDataFileItemRenderer();
    private MetaDataImageRenderer imageRenderer = new MetaDataImageRenderer();
    private MetaDataStringRenderer stringRenderer = new MetaDataStringRenderer();
    private MetaDataDoubleRenderer doubleRenderer = new MetaDataDoubleRenderer();
    private MetaDataIntegerRenderer numberRenderer = new MetaDataIntegerRenderer();
    //private TableImageItemEditor imageEditor = new TableImageItemEditor();
    private TableFileItemEditor metadataEditor = new TableFileItemEditor(new JTextField());
    private TableFileItemEditor fileEditor = new TableFileItemEditor(new JTextField());
    private TableRowSorter sorter;
    private RowEnableFilter rowEnableFilter = new RowEnableFilter();
    private JDialogColumnsSelector frameColumnsSelector = new JDialogColumnsSelector();

    /** Creates new form JFrameMetaData */
    public JFrameMetaData(String filename) {
        super();

        setTitle(getTitle(filename));

        initComponents();

        // Builds table.
        tableModel = new MetaDataTableModel(filename);

        // Build blocks list.
        String blocks[] = tableModel.getBlocks();
        for (int i = 0; i < blocks.length; i++) {
            jcbBlock.addItem(blocks[i]);
        }

        table = new JTable(tableModel) {

            //Implement table header tool tips.
            @Override
            protected JTableHeader createDefaultTableHeader() {
                return new JTableHeader(columnModel) {

                    @Override
                    public String getToolTipText(MouseEvent e) {
                        java.awt.Point p = e.getPoint();
                        int index = columnModel.getColumnIndexAtX(p.x);
                        int realIndex = columnModel.getColumn(index).getModelIndex();
                        return tableModel.getColumnName(realIndex);
                    }
                };
            }
        };

        rowEnableFilter = new RowEnableFilter();
        sorter = new TableRowSorter<MetaDataTableModel>(tableModel);
        sorter.setSortsOnUpdates(true);
        table.setRowSorter(sorter);

        table.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);

        table.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                tableMouseClicked(evt);
            }
        });
        table.addKeyListener(new KeyListener() {

            public void keyTyped(KeyEvent ke) {
            }

            public void keyPressed(KeyEvent ke) {
            }

            public void keyReleased(KeyEvent ke) {
                if (ke.getKeyCode() == KeyEvent.VK_DELETE) {
                    removeSelectedRows();
                }
            }
        });

        jsPanel.setViewportView(table);
//        setRowHeader();

        table.setColumnModel(columnModel);

        setRenderers();
        enableEditors(true);

        jcbBlock.setSelectedIndex(tableModel.getSelectedBlockIndex());
        jcbBlock.addItemListener(new java.awt.event.ItemListener() {

            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                if (evt.getStateChange() == ItemEvent.SELECTED) {
                    selectBlock(jcbBlock.getSelectedIndex());
                }
            }
        });

        pack();
        reloadTableData();
    }

    private String getTitle(String title) {
        int strlenght = title.length() / 4; // Approximated string graphical lenght.

        int toremove = strlenght - getWidth();

        int start = toremove > 0 ? toremove : 0;

        String newtitle = title.substring(start, title.length());
        String prefix = toremove > 0 ? "..." : "";

        return prefix + newtitle;
    }

    private void setRowHeader() {
        rowHeaderModel = new ImagesRowHeaderModel(table);

        rowHeader = new JList();
        rowHeader.setModel(rowHeaderModel);

        LookAndFeel.installColorsAndFont(rowHeader, "TableHeader.background",
                "TableHeader.foreground", "TableHeader.font");

//        rowHeader.setFixedCellHeight(imageRenderer.getCellHeight());

        rowHeader.setCellRenderer(rowHeaderRenderer);

        jsPanel.setRowHeaderView(rowHeader);

        rowHeader.addListSelectionListener(new ListSelectionListener() {

            public void valueChanged(ListSelectionEvent lse) {
                if (!lse.getValueIsAdjusting()) {
                    System.out.println("@TODO Selection rowHeader <-> table.");
                }
            }
        });
    }

    private void setRenderers() {
        table.setDefaultRenderer(GalleryImageItem.class, imageRenderer);
        table.setDefaultRenderer(TableMetaDataItem.class, fileRenderer);
        table.setDefaultRenderer(TableFileItem.class, fileRenderer);
        table.setDefaultRenderer(String.class, stringRenderer);
        table.setDefaultRenderer(Double.class, doubleRenderer);
        table.setDefaultRenderer(Integer.class, numberRenderer);
        //table.setDefaultRenderer(Boolean.class , booleanRenderer);
    }

    private void enableEditors(boolean enable) {
//        table.setDefaultEditor(GalleryImageItem.class, enable ? imageEditor : null);
        table.setDefaultEditor(TableMetaDataItem.class, enable ? metadataEditor : null);
        table.setDefaultEditor(TableFileItem.class, enable ? fileEditor : null);
    }

    private synchronized void updateTableStructure() {
        DEBUG.printMessage(" *** Updating table... " + System.currentTimeMillis());
        packColumns();
        packRows();

        table.getTableHeader().revalidate();
        table.revalidate();
        rowHeader.repaint();

//        jsPanel.setRowHeaderView(rowHeader);
    }

    private void packColumns() {
        for (int column = 0; column < table.getColumnCount(); column++) {
            int r = 0;

            // Get width of column header
            TableColumn tc = columnModel.getColumn(column);

            TableCellRenderer renderer = tc.getHeaderRenderer();
            if (renderer == null) {
                renderer = table.getTableHeader().getDefaultRenderer();
            }
            Component comp = renderer.getTableCellRendererComponent(
                    table, tc.getHeaderValue(), false, false, 0, 0);
            int headerWidth = comp.getPreferredSize().width;

            renderer = table.getCellRenderer(r, column);
            comp = renderer.getTableCellRendererComponent(
                    table, table.getValueAt(r, column), false, false, r, column);
            int width = comp.getPreferredSize().width;

            tc.setPreferredWidth(Math.max(headerWidth, width));
        }
    }

    // The height of each row is set to the preferred height of the tallest cell in that row.
    private void packRows() {
//        int cellHeight = MetaDataImageRenderer.DEFAULT_CELL_HEIGHT;

        int row = 0;
        int height = -1;
        //table.setRowHeight(row, cellHeight);

        for (int column = 0; column < tableModel.getColumnCount(); column++) {
            TableCellRenderer renderer = table.getCellRenderer(row, column);
            Component comp = renderer.getTableCellRendererComponent(
                    table, table.getValueAt(row, column), false, false, row, column);
            height = Math.max(height, comp.getPreferredSize().height);
        }
//            TableColumn tcolumn = columnModel.getColumn(column);

        // Sets width
        table.setRowHeight(height);

        //  Row header.
        rowHeader.setFixedCellHeight(height);
    }

    private void tableMouseClicked(java.awt.event.MouseEvent evt) {
        int view_row = table.rowAtPoint(evt.getPoint());
        int view_col = table.columnAtPoint(evt.getPoint());

        if (view_row >= 0 && view_col >= 0) {
            if (SwingUtilities.isLeftMouseButton(evt)) {
                //if (jcbRenderImages.isSelected()) { // Opens images just when rendering. Otherwise, edits it.
                if (evt.getClickCount() > 1) {
                    Object item = table.getValueAt(view_row, view_col);

                    if (item instanceof GalleryImageItem) {
                        openXmippImage((GalleryImageItem) item);
                    } else if (item instanceof TableMetaDataItem) {
                        openMetaDataFile((TableMetaDataItem) item);
                    } else if (item instanceof TableFileItem) {
                        openXmippTextFile((TableFileItem) item);
                    }
                }
                //}
            } else if (SwingUtilities.isRightMouseButton(evt)) {
                table.setRowSelectionInterval(view_row, view_row);
                table.setColumnSelectionInterval(view_col, view_col);
            }
        }
        updateTableStructure();
    }

    private void openXmippImage(GalleryImageItem item) {
        ImagesWindowFactory.captureFrame(item.getImagePlus());
    }

    private void openXmippTextFile(TableFileItem item) {
        JDialogTextFile frame = new JDialogTextFile(item.getPath(), item.getOriginalValue());
        frame.setLocationRelativeTo(this);
        frame.setVisible(true);
    }

    private void openMetaDataFile(TableMetaDataItem item) {
        JFrameMetaData frameMetaData = new JFrameMetaData(item.getPath());
        setTitle(getTitle(item.getOriginalValue()));
        frameMetaData.setLocationRelativeTo(this);
        frameMetaData.setVisible(true);
    }

    private void save(String fileName) {
        if (tableModel.save(fileName)) {
            IJ.showMessage(Labels.MESSAGE_FILE_SAVED + fileName);
        } else {
            IJ.error("File NOT saved: " + fileName);
        }
    }

    private void enableAll(boolean enable) {
        tableModel.enableAllRows(enable);

        updateTableStructure();
    }

    private void hideDisabled(boolean hide) {
        sorter.setRowFilter(hide ? rowEnableFilter : null);

        updateTableStructure();
    }

    private void enableImagesRendering(boolean renderImages) {
        imageRenderer.setRenderImages(renderImages);
        enableEditors(!renderImages);

        updateTableStructure();
    }

    private void reloadTableData() {
        DEBUG.printMessage(" *** Reloading table: " + System.currentTimeMillis());

        columnModel.clear();
        tableModel.reload();

        canRender = tableModel.containsImages();
        if (!canRender) {
            jcbRenderImages.setSelected(false);
        }
        jcbRenderImages.setEnabled(canRender);
        jbSend2Gallery.setEnabled(canRender);

        setRowHeader();

        updateTableStructure();
    }

    private void selectBlock(int block) {
        tableModel.selectBlock(block);
        jbSend2Gallery.setEnabled(tableModel.containsImages());
        reloadTableData();
    }

    private void setColumnsVisible(boolean states[]) {
        for (int i = 0; i < columnModel.getColumnCount(false); i++) {
            columnModel.setColumnVisible(columnModel.getColumnByModelIndex(i), states[i]);
        }
    }

    private boolean[] getColumnsStates() {
        boolean states[] = new boolean[columnModel.getColumnCount(false)];

        for (int i = 0; i < states.length; i++) {
            states[i] = columnModel.isColumnVisible(columnModel.getColumnByModelIndex(i));
        }

        return states;
    }

    private void removeSelectedRows() {
        int i;

        while ((i = table.getSelectedRow()) != -1) {
            tableModel.removeRow(table.convertRowIndexToModel(i));
        }

        updateTableStructure();
    }

    public void setRenderImages(boolean render_images) {
        if (canRender) {
            jcbRenderImages.setSelected(render_images);
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

        toolBar = new javax.swing.JToolBar();
        jbSave = new javax.swing.JButton();
        jbReload = new javax.swing.JButton();
        jbSend2Gallery = new javax.swing.JButton();
        jpCenter = new javax.swing.JPanel();
        jpControls = new javax.swing.JPanel();
        jlBlock = new javax.swing.JLabel();
        jcbBlock = new javax.swing.JComboBox();
        jcbEnableAll = new javax.swing.JCheckBox();
        jcbHideDisabled = new javax.swing.JCheckBox();
        jcbRenderImages = new javax.swing.JCheckBox();
        bHideColumns = new javax.swing.JButton();
        jlStatus = new javax.swing.JLabel();
        jsPanel = new javax.swing.JScrollPane();

        toolBar.setRollover(true);

        jbSave.setText(METADATA_LABELS.SAVE);
        jbSave.setFocusable(false);
        jbSave.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbSave.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbSave.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbSaveActionPerformed(evt);
            }
        });
        toolBar.add(jbSave);

        jbReload.setText(METADATA_LABELS.RELOAD);
        jbReload.setFocusable(false);
        jbReload.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbReload.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbReload.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbReloadActionPerformed(evt);
            }
        });
        toolBar.add(jbReload);

        jbSend2Gallery.setText(METADATA_LABELS.OPEN_AS_GALLERY);
        jbSend2Gallery.setFocusable(false);
        jbSend2Gallery.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbSend2Gallery.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbSend2Gallery.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbSend2GalleryActionPerformed(evt);
            }
        });
        toolBar.add(jbSend2Gallery);

        getContentPane().add(toolBar, java.awt.BorderLayout.PAGE_START);

        jpCenter.setLayout(new java.awt.BorderLayout());

        jpControls.setLayout(new javax.swing.BoxLayout(jpControls, javax.swing.BoxLayout.LINE_AXIS));

        jlBlock.setText(METADATA_LABELS.BLOCK);
        jpControls.add(jlBlock);
        jpControls.add(jcbBlock);

        jcbEnableAll.setText(METADATA_LABELS.ENABLE_ALL);
        jcbEnableAll.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jcbEnableAllItemStateChanged(evt);
            }
        });
        jpControls.add(jcbEnableAll);

        jcbHideDisabled.setText(METADATA_LABELS.HIDE_DISABLED);
        jcbHideDisabled.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jcbHideDisabledItemStateChanged(evt);
            }
        });
        jpControls.add(jcbHideDisabled);

        jcbRenderImages.setText(METADATA_LABELS.RENDER_IMAGES);
        jcbRenderImages.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jcbRenderImagesItemStateChanged(evt);
            }
        });
        jpControls.add(jcbRenderImages);

        bHideColumns.setText(METADATA_LABELS.SELECT_COLUMNS);
        bHideColumns.setFocusable(false);
        bHideColumns.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        bHideColumns.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        bHideColumns.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                bHideColumnsActionPerformed(evt);
            }
        });
        jpControls.add(bHideColumns);
        jpControls.add(jlStatus);

        jpCenter.add(jpControls, java.awt.BorderLayout.PAGE_START);
        jpCenter.add(jsPanel, java.awt.BorderLayout.CENTER);

        getContentPane().add(jpCenter, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jbSaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbSaveActionPerformed
        if (fc == null) {
            // Sets path and filename automatically.
            fc = new JFileChooser();
            fc.setSelectedFile(new File(tableModel.getFilename()));
        }

        if (fc.showSaveDialog(this) != JFileChooser.CANCEL_OPTION) {
            boolean response = true;
            if (fc.getSelectedFile().exists()) {
                response = JOptionPane.showConfirmDialog(null,
                        "Overwrite existing file?", "Confirm Overwrite",
                        JOptionPane.OK_CANCEL_OPTION,
                        JOptionPane.QUESTION_MESSAGE) == JOptionPane.OK_OPTION;
            }

            if (response) {
                save(fc.getSelectedFile().getAbsolutePath());
            }
        }
    }//GEN-LAST:event_jbSaveActionPerformed

    private void jcbEnableAllItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jcbEnableAllItemStateChanged
        enableAll(jcbEnableAll.isSelected());
    }//GEN-LAST:event_jcbEnableAllItemStateChanged

    private void jcbHideDisabledItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jcbHideDisabledItemStateChanged
        hideDisabled(jcbHideDisabled.isSelected());
}//GEN-LAST:event_jcbHideDisabledItemStateChanged

    private void jbReloadActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbReloadActionPerformed
        reloadTableData();
}//GEN-LAST:event_jbReloadActionPerformed

    private void jcbRenderImagesItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jcbRenderImagesItemStateChanged
        enableImagesRendering(jcbRenderImages.isSelected());
    }//GEN-LAST:event_jcbRenderImagesItemStateChanged

    private void bHideColumnsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_bHideColumnsActionPerformed
        frameColumnsSelector.setModel(tableModel);
        frameColumnsSelector.setStates(getColumnsStates());
        frameColumnsSelector.setLocationRelativeTo(this);
        frameColumnsSelector.setVisible(true);

        if (frameColumnsSelector.hasSelectedOK()) {
            setColumnsVisible(frameColumnsSelector.getStates());
        }
    }//GEN-LAST:event_bHideColumnsActionPerformed

private void jbSend2GalleryActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbSend2GalleryActionPerformed
    ImagesWindowFactory.openFileAsGallery(tableModel.getPath());
}//GEN-LAST:event_jbSend2GalleryActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton bHideColumns;
    private javax.swing.JButton jbReload;
    private javax.swing.JButton jbSave;
    private javax.swing.JButton jbSend2Gallery;
    private javax.swing.JComboBox jcbBlock;
    private javax.swing.JCheckBox jcbEnableAll;
    private javax.swing.JCheckBox jcbHideDisabled;
    private javax.swing.JCheckBox jcbRenderImages;
    private javax.swing.JLabel jlBlock;
    private javax.swing.JLabel jlStatus;
    private javax.swing.JPanel jpCenter;
    private javax.swing.JPanel jpControls;
    private javax.swing.JScrollPane jsPanel;
    private javax.swing.JToolBar toolBar;
    // End of variables declaration//GEN-END:variables
}
