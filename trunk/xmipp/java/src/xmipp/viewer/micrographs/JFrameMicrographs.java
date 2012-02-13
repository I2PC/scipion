/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * JFrameMicrographs.java
 *
 * Created on Oct 25, 2010, 10:50:01 AM
 */
package xmipp.viewer.micrographs;

import xmipp.utils.DEBUG;
import xmipp.utils.Resources;
import xmipp.utils.Labels;
import xmipp.viewer.gallery.models.GalleryRowHeaderModel;
import xmipp.viewer.gallery.renderers.RowHeaderRenderer;
import xmipp.viewer.imageitems.tableitems.GalleryImageItem;
import xmipp.viewer.windows.ImagesWindowFactory;
import ij.IJ;
import ij.ImagePlus;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.util.ArrayList;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.JSeparator;
import javax.swing.JTable;
import javax.swing.LookAndFeel;
import javax.swing.SwingUtilities;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import javax.swing.table.TableRowSorter;
import xmipp.viewer.metadata.images.TableFileItem;
import xmipp.viewer.metadata.images.TableMetaDataItem;
import xmipp.viewer.metadata.models.XTableColumnModel;
import xmipp.viewer.metadata.renderers.MetaDataDoubleRenderer;
import xmipp.viewer.metadata.renderers.MetaDataFileItemRenderer;
import xmipp.viewer.metadata.renderers.MetaDataImageRenderer;
import xmipp.viewer.metadata.renderers.MetaDataIntegerRenderer;
import xmipp.viewer.metadata.renderers.MetaDataStringRenderer;
import xmipp.viewer.micrographs.ctf.tasks.TasksEngine;
import xmipp.viewer.micrographs.ctf.tasks.iCTFGUI;
import xmipp.viewer.micrographs.filters.EnableFilter;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameMicrographs extends JFrame implements iCTFGUI {

    private final static int CTF_IMAGE_COLUMN = 3;
    private JTable table;
    private MicrographsTableModel tableModel;
    private GalleryRowHeaderModel rowHeaderModel;
    private XTableColumnModel columnModel = new XTableColumnModel();
    private JPopUpMenuMicrograph jPopupMenu = new JPopUpMenuMicrograph();
    private JFileChooser fc = new JFileChooser();
    private JList rowHeader;
    private MetaDataFileItemRenderer fileRenderer = new MetaDataFileItemRenderer();
    private MetaDataImageRenderer imageRenderer = new MetaDataImageRenderer();
    private MetaDataStringRenderer stringRenderer = new MetaDataStringRenderer();
    private MetaDataDoubleRenderer doubleRenderer = new MetaDataDoubleRenderer();
    private MetaDataIntegerRenderer numberRenderer = new MetaDataIntegerRenderer();
    private TableRowSorter sorter;
    private EnableFilter enableFilter;
    private TasksEngine tasksEngine = new TasksEngine(this);

    /** Creates new form JFrameMicrographs */
    public JFrameMicrographs(String filename) {
        super();

        setTitle(filename);//ImagesWindowFactory.getSortTitle(filename, getWidth(),
        //getGraphics().getFontMetrics()));

        initComponents();

        // Builds table.
        tableModel = new MicrographsTableModel(filename);

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

        // Sets sorter and filter.
/*        sorter = new TableRowSorter<MicrographsTableModel>(tableModel);
        sorter.setSortsOnUpdates(true);
        enableFilter = new EnableFilter();
        sorter.setRowFilter(enableFilter);*/

        enableFilter = new EnableFilter();
        setNewTableRowSorter();
        table.setRowSorter(sorter);

        table.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);

        table.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                tableMouseClicked(evt);
            }
        });

        jsPanel.setViewportView(table);

        table.setColumnModel(columnModel);
        table.createDefaultColumnsFromModel();

        setRenderers();

        // places Defocus columns next to images.
        table.moveColumn(table.getColumnCount() - 2, MicrographsTableModel.INDEX_DEFOCUS_U_COL);
        table.moveColumn(table.getColumnCount() - 1, MicrographsTableModel.INDEX_DEFOCUS_V_COL);

        hideColumns();
        setRowHeader();

        updateTableStructure();
        autoSortTable(MicrographsTableModel.INDEX_AUTOSORT_COLUMN);

        pack();
    }

    private void autoSortTable(int column) {
        int view_column = table.convertColumnIndexToView(column);

        if (column < table.getColumnCount()) {
            table.getRowSorter().toggleSortOrder(view_column);
        }
    }

    private void setNewTableRowSorter() {
        sorter = new TableRowSorter<MicrographsTableModel>(tableModel);
        sorter.setSortsOnUpdates(true);
        sorter.setRowFilter(enableFilter);
        table.setRowSorter(sorter);
    }

    public String getFilename() {
        return tableModel.getFilename();
    }

    public void setRowBusy(int row) {
        setRunning(true);
        tableModel.setRowBusy(row);
    }

    public void setRowIdle(int row) {
        tableModel.setRowIdle(row);

        jPopupMenu.refresh();   // If menu is showing it will be refreshed.
    }

    private void hideColumns() {
        ArrayList<Integer> columns = tableModel.getColumnsToHide();

        for (int i = 0; i < columns.size(); i++) {
            TableColumn column = columnModel.getColumnByModelIndex(columns.get(i));
            columnModel.setColumnVisible(column, false);
        }
    }

    private void setRowHeader() {
        rowHeaderModel = new GalleryRowHeaderModel(tableModel.getRowCount());

        rowHeader = new JList();
        rowHeader.setModel(rowHeaderModel);

        LookAndFeel.installColorsAndFont(rowHeader, "TableHeader.background",
                "TableHeader.foreground", "TableHeader.font");

        rowHeader.setCellRenderer(new RowHeaderRenderer());

        jsPanel.setRowHeaderView(rowHeader);
    }

    private void setRenderers() {
        table.setDefaultRenderer(GalleryImageItem.class, imageRenderer);
        table.setDefaultRenderer(TableMetaDataItem.class, fileRenderer);
        table.setDefaultRenderer(TableFileItem.class, fileRenderer);
        table.setDefaultRenderer(String.class, stringRenderer);
        table.setDefaultRenderer(Double.class, doubleRenderer);
        table.setDefaultRenderer(Integer.class, numberRenderer);
    }

    // Sets the preferred width of the visible column specified by vColIndex. The column
    // will be just wide enough to show the column head and the widest cell in the column.
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
        int row = 0;
        int height = 0;
        //table.setRowHeight(row, cellHeight);

        for (int column = 0; column < table.getColumnCount(); column++) {

            TableCellRenderer renderer = table.getCellRenderer(row, column);
            Component comp = renderer.getTableCellRendererComponent(
                    table, table.getValueAt(row, column), false, false, row, column);

            height = Math.max(height, comp.getPreferredSize().height);
        }

        // Sets width
        table.setRowHeight(height);

        // Row header.
        rowHeader.setFixedCellHeight(height);
    }

    private void tableMouseClicked(java.awt.event.MouseEvent evt) {
        int view_row = table.rowAtPoint(evt.getPoint());
        int view_col = table.columnAtPoint(evt.getPoint());
        int model_row = table.convertRowIndexToModel(view_row);
        int model_col = table.convertColumnIndexToModel(view_col);

        if (SwingUtilities.isLeftMouseButton(evt)) {
            if (evt.getClickCount() > 1) {
                Object item = table.getValueAt(view_row, view_col);

                if (item instanceof GalleryImageItem) {
                    GalleryImageItem tableItem = (GalleryImageItem) item;
                    ImagePlus ip = tableItem.getImagePlus();

                    if (ip != null) {
                        ip.setTitle(tableItem.getOriginalValue());

                        ImagesWindowFactory.captureFrame(ip);
                    }
                }
            } else {    // Single click
                if (model_col == MicrographsTableModel.INDEX_ENABLED) {
                    int min = table.getSelectionModel().getMinSelectionIndex();
                    int max = table.getSelectionModel().getMaxSelectionIndex();

                    // Multiple rows selection.
                    if (evt.isShiftDown()) {
                        // Last selection value will be the value for all the selected items.
                        Boolean value = (Boolean) tableModel.getValueAt(model_row, model_col);

                        for (int i = min; i <= max; i++) {
                            tableModel.setValueAt(value, table.convertRowIndexToModel(i),
                                    model_col);
                        }
                    }

                    updateTableStructure();
                }
            }
        } else if (SwingUtilities.isRightMouseButton(evt)) {
            table.setRowSelectionInterval(view_row, view_row);
            table.setColumnSelectionInterval(view_col, view_col);

            jPopupMenu.setCell(model_row, model_col);

            jPopupMenu.refresh();
            jPopupMenu.show(evt.getComponent(), evt.getX(), evt.getY());
        }
    }

    private void save(String fileName) {
        if (tableModel.save(fileName)) {
            IJ.showMessage(Labels.MESSAGE_FILE_SAVED + fileName);
        }
    }

    private void enableAll(boolean enable) {
        for (int i = 0; i < tableModel.getRowCount(); i++) {
            tableModel.setValueAt(enable, i, MicrographsTableModel.ENABLED_COLUMN_INDEX);
        }

        tableModel.fireTableRowsUpdated(0, tableModel.getRowCount() - 1);

        updateTableStructure();
    }

    private void setFiltering(boolean filtering) {
        enableFilter.setFiltering(filtering);

        table.setColumnModel(columnModel);  // Re sets column model to hide columns.

        tableModel.fireTableDataChanged();

        updateTableStructure();
    }

    private void showCTFImage(GalleryImageItem item, String CTFFilename,
            String PSDfilename, String MicrographFilename, int row) {
        ImagesWindowFactory.openCTFImage(item.getImagePlus(), CTFFilename,
                PSDfilename, tasksEngine, MicrographFilename, row);
    }

    public void setRunning(boolean running) {
        jlStatus.setIcon(running ? Resources.WAIT_ICON : null);
    }

    public void done() {
        setRunning(false);

        refreshTableData();
    }

    private void refreshTableData() {
        DEBUG.printMessage(" *** Refreshing table [AUTO]: " + System.currentTimeMillis());

        boolean enabled[] = tableModel.getEnabledRows();   // Gets currently enabled rows...
        tableModel.reload();
        tableModel.setEnabledRows(enabled);    // ...sets previously enabled.
        table.setColumnModel(columnModel);  // Re sets column model to hide columns.

        updateTableStructure();
    }

    private void reloadTableData() {
        DEBUG.printMessage(" *** Refreshing table [RELOAD]: " + System.currentTimeMillis());

        tableModel.reload();
        table.setColumnModel(columnModel);  // Re sets column model to hide columns.
        setNewTableRowSorter();

        updateTableStructure();
    }

    private synchronized void updateTableStructure() {
        ImagesWindowFactory.blockGUI(getRootPane(), "Updating table...");

        Thread t = new Thread(new Runnable() {

            public void run() {
                packColumns();
                packRows();

                rowHeader.repaint();
                table.getTableHeader().repaint();

                ImagesWindowFactory.releaseGUI(getRootPane());
            }
        });

        t.start();
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
        bSave = new javax.swing.JButton();
        bReload = new javax.swing.JButton();
        jpCenter = new javax.swing.JPanel();
        jpCheckAll = new javax.swing.JPanel();
        jcbEnableAll = new javax.swing.JCheckBox();
        jcbFilterEnabled = new javax.swing.JCheckBox();
        jlStatus = new javax.swing.JLabel();
        jsPanel = new javax.swing.JScrollPane();

        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowOpened(java.awt.event.WindowEvent evt) {
                formWindowOpened(evt);
            }
        });

        toolBar.setRollover(true);

        bSave.setText(Labels.BUTTON_SAVE);
        bSave.setFocusable(false);
        bSave.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        bSave.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        bSave.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                bSaveActionPerformed(evt);
            }
        });
        toolBar.add(bSave);

        bReload.setText(Labels.BUTTON_RELOAD_GALLERY);
        bReload.setFocusable(false);
        bReload.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        bReload.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        bReload.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                bReloadActionPerformed(evt);
            }
        });
        toolBar.add(bReload);

        getContentPane().add(toolBar, java.awt.BorderLayout.PAGE_START);

        jpCenter.setLayout(new java.awt.BorderLayout());

        jpCheckAll.setLayout(new javax.swing.BoxLayout(jpCheckAll, javax.swing.BoxLayout.LINE_AXIS));

        jcbEnableAll.setText(Labels.LABEL_GALLERY_ENABLE_ALL);
        jcbEnableAll.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jcbEnableAllItemStateChanged(evt);
            }
        });
        jpCheckAll.add(jcbEnableAll);

        jcbFilterEnabled.setText(Labels.LABEL_GALLERY_HIDE_DISABLED);
        jcbFilterEnabled.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jcbFilterEnabledItemStateChanged(evt);
            }
        });
        jpCheckAll.add(jcbFilterEnabled);
        jpCheckAll.add(jlStatus);

        jpCenter.add(jpCheckAll, java.awt.BorderLayout.PAGE_START);
        jpCenter.add(jsPanel, java.awt.BorderLayout.CENTER);

        getContentPane().add(jpCenter, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void bSaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_bSaveActionPerformed
        // Sets path and filename automatically.
        fc.setSelectedFile(new File(tableModel.getFilename()));

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
    }//GEN-LAST:event_bSaveActionPerformed

    private void jcbEnableAllItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jcbEnableAllItemStateChanged
        enableAll(jcbEnableAll.isSelected());
    }//GEN-LAST:event_jcbEnableAllItemStateChanged

    private void jcbFilterEnabledItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jcbFilterEnabledItemStateChanged
        setFiltering(jcbFilterEnabled.isSelected());
}//GEN-LAST:event_jcbFilterEnabledItemStateChanged

    private void bReloadActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_bReloadActionPerformed
        reloadTableData();
}//GEN-LAST:event_bReloadActionPerformed

private void formWindowOpened(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowOpened
    ImagesWindowFactory.setConvenientSize(this);
}//GEN-LAST:event_formWindowOpened
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton bReload;
    private javax.swing.JButton bSave;
    private javax.swing.JCheckBox jcbEnableAll;
    private javax.swing.JCheckBox jcbFilterEnabled;
    private javax.swing.JLabel jlStatus;
    private javax.swing.JPanel jpCenter;
    private javax.swing.JPanel jpCheckAll;
    private javax.swing.JScrollPane jsPanel;
    private javax.swing.JToolBar toolBar;
    // End of variables declaration//GEN-END:variables

    class JPopUpMenuMicrograph extends JPopupMenu {

        private int row, col;   // Current cell where menu is displayed.
        private JMenuItem jmiShowCTF = new JMenuItem(Labels.LABEL_SHOW_CTF_INFO);
        private JMenuItem jmiExtractColumnEnabled = new JMenuItem(Labels.LABEL_EXTRACT_COLUMN_ENABLED);
        private JMenuItem jmiExtractColumnAll = new JMenuItem(Labels.LABEL_EXTRACT_COLUMN_ALL);
        private JMenuItem jmiRecalculateCTF = new JMenuItem(Labels.LABEL_RECALCULATE_CTF);
        private JMenuItem jmiViewCTFProfile = new JMenuItem(Labels.LABEL_VIEW_CTF_PROFILE);

        public JPopUpMenuMicrograph() {
            super();

            add(jmiShowCTF);
            add(new JSeparator());
            add(jmiExtractColumnEnabled);
            add(jmiExtractColumnAll);
            add(new JSeparator());
            add(jmiViewCTFProfile);
            add(jmiRecalculateCTF);

            jmiShowCTF.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    showCTFFile(table.getSelectedRow());
                }
            });

            jmiExtractColumnEnabled.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    extractColumn(table.getSelectedColumn(), true);
                }
            });

            jmiExtractColumnAll.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    extractColumn(table.getSelectedColumn(), false);
                }
            });

            jmiViewCTFProfile.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    showCTFProfile();
                }
            });

            jmiRecalculateCTF.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    showRecalculateCTFWindow();
                }
            });
        }

        private void setCell(int row, int col) {
            this.row = row;
            this.col = col;
        }

        private void refresh() {
            boolean hasCtfData = tableModel.hasCtfData();
            jmiShowCTF.setEnabled(hasCtfData);
            jmiViewCTFProfile.setEnabled(hasCtfData);

            // Column extraction only allowed for images
            jmiExtractColumnEnabled.setEnabled(tableModel.getValueAt(row, col) instanceof GalleryImageItem);
            jmiExtractColumnAll.setEnabled(tableModel.getValueAt(row, col) instanceof GalleryImageItem);

            boolean busy = tableModel.isRowBusy(row);
            jmiRecalculateCTF.setIcon(busy ? Resources.WAIT_MENU_ICON : null);
            jmiRecalculateCTF.setEnabled(!busy && hasCtfData);

//            repaint();
            updateUI();

            pack();
        }

        private void showCTFFile(int row) {
            String CTFfile = tableModel.getCTFfile(row);

            ImagesWindowFactory.openFileAsText(CTFfile, this);
        }

        private void extractColumn(int column, boolean onlyenabled) {
            if (onlyenabled) {
                ImagesWindowFactory.openFilesAsGallery(tableModel.extractColumn(column, onlyenabled), true);
            } else {
                ImagesWindowFactory.openGallery(
                        tableModel.extractColumn(column, onlyenabled),
                        tableModel.getEnabledRows());
            }
        }

        private void showCTFProfile() {
            String CTFFilename = tableModel.getCTFfile(row);
            String displayFilename = tableModel.getCTFDisplayfile(row);
            String PSDFilename = tableModel.getPSDfile(row);

            ImagePlus ip = IJ.openImage(displayFilename);

            ImagesWindowFactory.openCTFWindow(ip, CTFFilename, PSDFilename);
        }

        private void showRecalculateCTFWindow() {
            Object item = tableModel.getValueAt(row, CTF_IMAGE_COLUMN);

            if (item instanceof GalleryImageItem) {
                showCTFImage((GalleryImageItem) item, tableModel.getCTFfile(row),
                        tableModel.getPSDfile(row), tableModel.getFilename(), row);
            }
        }
    }
}
