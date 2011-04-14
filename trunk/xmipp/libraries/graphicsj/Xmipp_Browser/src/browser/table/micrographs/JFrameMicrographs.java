/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * JFrameMicrographs.java
 *
 * Created on Oct 25, 2010, 10:50:01 AM
 */
package browser.table.micrographs;

import browser.table.micrographs.ctf.JFrameCTF;
import browser.LABELS;
import browser.imageitems.TableImageItem;
import browser.table.ImagesRowHeaderModel;
import browser.table.micrographs.filters.EnableFilter;
import browser.table.micrographs.renderers.MicrographDoubleRenderer;
import browser.table.micrographs.renderers.MicrographFileNameRenderer;
import browser.table.micrographs.renderers.MicrographImageRenderer;
import browser.table.renderers.RowHeaderRenderer;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.LinkedList;
import java.util.Vector;
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
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import javax.swing.table.TableRowSorter;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameMicrographs extends JFrame implements iMicrographsGUI {

    private final static int CTF_IMAGE_COLUMN = 3;
    private JTable table;
    private MicrographsTableModel tableModel;
    private ImagesRowHeaderModel rowHeaderModel;
    private XTableColumnModel columnModel = new XTableColumnModel();
    private JPopUpMenuMicrograph jPopupMenu = new JPopUpMenuMicrograph();
    private JFileChooser fc = new JFileChooser();
    private JList rowHeader;
    private MicrographFileNameRenderer fileNameRenderer = new MicrographFileNameRenderer();
    private MicrographImageRenderer imageRenderer = new MicrographImageRenderer();
    private MicrographDoubleRenderer doubleRenderer = new MicrographDoubleRenderer();
    private JFrameCTF frameCTF = new JFrameCTF();
    private TableRowSorter sorter;
    private EnableFilter enableFilter;

    /** Creates new form JFrameMicrographs */
    public JFrameMicrographs(String filename) {
        super();

        setTitle(getTitle(filename));

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
        sorter = new TableRowSorter<MicrographsTableModel>(tableModel);
        enableFilter = new EnableFilter();
        sorter.setRowFilter(enableFilter);
        table.setRowSorter(sorter);

        //table.setAutoCreateRowSorter(true);
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

        updateTable();
        autoSortTable(MicrographsTableModel.INDEX_COMBINED_COLUMN);

        pack();
    }

    private void autoSortTable(int column) {
        int view_column = table.convertColumnIndexToView(column);

        if (column < table.getColumnCount()) {
            table.getRowSorter().toggleSortOrder(view_column);
        }
    }

    private String getTitle(String title) {
        int strlenght = title.length() / 4; // Approximated string graphical lenght.

        int toremove = strlenght - getWidth();

        int start = toremove > 0 ? toremove : 0;

        String newtitle = title.substring(start, title.length());
        String prefix = toremove > 0 ? "..." : "";

        return prefix + newtitle;
    }

    private void hideColumns() {
        Vector<Integer> columns = tableModel.getColumnsToHide();

        for (int i = 0; i < columns.size(); i++) {
            TableColumn column = columnModel.getColumnByModelIndex(columns.elementAt(i));
            columnModel.setColumnVisible(column, false);
        }
    }

    private void setRowHeader() {
        rowHeaderModel = new ImagesRowHeaderModel(table);

        rowHeader = new JList();
        rowHeader.setModel(rowHeaderModel);

        LookAndFeel.installColorsAndFont(rowHeader, "TableHeader.background",
                "TableHeader.foreground", "TableHeader.font");

        rowHeader.setFixedCellHeight(MicrographImageRenderer.CELL_HEIGHT);

        rowHeader.setCellRenderer(new RowHeaderRenderer());

        jsPanel.setRowHeaderView(rowHeader);
    }

    private void setRenderers() {
        table.setDefaultRenderer(TableImageItem.class, imageRenderer);
        table.setDefaultRenderer(Double.class, doubleRenderer);

        // Image is quite big, so don't show it by default.
        table.getColumnModel().getColumn(MicrographsTableModel.INDEX_IMAGE).
                setCellRenderer(fileNameRenderer);
    }

    private void updateTable() {
        packColumns();
        packRows();
        rowHeader.repaint();
    }

    private void packColumns() {
        for (int i = 0; i < table.getColumnCount(); i++) {
            packColumn(table, i);
        }
    }

    // Sets the preferred width of the visible column specified by vColIndex. The column
    // will be just wide enough to show the column head and the widest cell in the column.
    private void packColumn(JTable table, int column) {
        DefaultTableColumnModel colModel = (DefaultTableColumnModel) table.getColumnModel();
        TableColumn col = colModel.getColumn(column);
        int width = 0;

        // Get width of column header
        TableCellRenderer renderer = col.getHeaderRenderer();
        if (renderer == null) {
            renderer = table.getTableHeader().getDefaultRenderer();
        }

        width = MicrographImageRenderer.CELL_WIDTH_MIN;

        // Get maximum width of column data
        for (int r = 0; r < table.getRowCount(); r++) {
            renderer = table.getCellRenderer(r, column);
            Component comp = renderer.getTableCellRendererComponent(
                    table, table.getValueAt(r, column), false, false, r, column);
            width = Math.max(width, comp.getPreferredSize().width);
        }

        // Set the width
        col.setPreferredWidth(width);
    }

    // The height of each row is set to the preferred height of the tallest cell in that row.
    private void packRows() {
        for (int r = 0; r < table.getRowCount(); r++) {
            packRows(table, r);
        }
    }

    // For each row >= start and < end, the height of a row is set to
    // the preferred height of the tallest cell in that row.
    private void packRows(JTable table, int row) {
        // Now set the row height using the preferred height
        if (table.getRowHeight(row) != MicrographImageRenderer.CELL_HEIGHT) {
            table.setRowHeight(row, MicrographImageRenderer.CELL_HEIGHT);
        }
    }

    private void tableMouseClicked(java.awt.event.MouseEvent evt) {
        int view_row = table.rowAtPoint(evt.getPoint());
        int view_col = table.columnAtPoint(evt.getPoint());

        if (SwingUtilities.isLeftMouseButton(evt)) {
            if (evt.getClickCount() > 1) {
                Object item = table.getValueAt(view_row, view_col);//getModel().getValueAt(row, col);

                if (item instanceof TableImageItem) {
                    ImagesWindowFactory.openImage(((TableImageItem) item).getImagePlus());
                }
            } else {    // Single click
                int min = table.getSelectionModel().getMinSelectionIndex();
                int max = table.getSelectionModel().getMaxSelectionIndex();

                if (evt.isShiftDown()) {
                    // Last selection value will be the value for all the selected items.
                    Boolean value = (Boolean) table.getValueAt(view_row, 0);//getModel().getValueAt(row, 0);

                    for (int i = min; i <= max; i++) {
                        table.setValueAt(value, i, 0);
                    }
                }

                if (enableFilter.isFiltering()) {
                    tableModel.fireTableDataChanged();
                } else {
                    tableModel.fireTableRowsUpdated(min, max);
                }
                updateTable();
            }
        } else if (SwingUtilities.isRightMouseButton(evt)) {
            table.setRowSelectionInterval(view_row, view_row);
            table.setColumnSelectionInterval(view_col, view_col);

            // Column extraction only allowed for images
            jPopupMenu.refreshItems(view_row, view_col);
            jPopupMenu.show(evt.getComponent(), evt.getX(), evt.getY());
        }
    }

    private void save(String fileName) {
        if (tableModel.save(fileName)) {
            IJ.showMessage("File saved: " + fileName);
        }
    }

    private void enableAll(boolean enable) {
        for (int i = 0; i < tableModel.getRowCount(); i++) {
            tableModel.setValueAt(enable, i, MicrographsTableModel.ENABLED_COLUMN_INDEX);
        }

        if (enableFilter.isFiltering()) {
            tableModel.fireTableDataChanged();
        } else {
            tableModel.fireTableRowsUpdated(0, tableModel.getRowCount() - 1);
        }
        updateTable();
    }

    private void setFiltering(boolean filtering) {
        enableFilter.setFiltering(filtering);

        tableModel.fireTableDataChanged();
        updateTable();
    }

    private void showCTFImage(TableImageItem item, String CTFFilename) {
        ImagesWindowFactory.openCTFImage(item.getImagePlus(), CTFFilename, this);
    }

    public void refresh() {
        System.out.println("@TODO RELOAD TABLE!");
/*        tableModel.reload();

        hideColumns();
        setRowHeader();

        updateTable();
        autoSortTable(MicrographsTableModel.INDEX_COMBINED_COLUMN);
*/
        // = new MicrographsTableModel(tableModel.getMicrographFilename());
//        table.setModel(tableModel);
//        updateTable();
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
        jpCenter = new javax.swing.JPanel();
        jpCheckAll = new javax.swing.JPanel();
        jcbEnableAll = new javax.swing.JCheckBox();
        jcbFilterEnabled = new javax.swing.JCheckBox();
        jsPanel = new javax.swing.JScrollPane();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        toolBar.setRollover(true);

        bSave.setText(LABELS.BUTTON_SAVE);
        bSave.setFocusable(false);
        bSave.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        bSave.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        bSave.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                bSaveActionPerformed(evt);
            }
        });
        toolBar.add(bSave);

        getContentPane().add(toolBar, java.awt.BorderLayout.PAGE_START);

        jpCenter.setLayout(new java.awt.BorderLayout());

        jpCheckAll.setLayout(new javax.swing.BoxLayout(jpCheckAll, javax.swing.BoxLayout.LINE_AXIS));

        jcbEnableAll.setText(LABELS.LABEL_TABLE_ENABLE_ALL);
        jcbEnableAll.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jcbEnableAllItemStateChanged(evt);
            }
        });
        jpCheckAll.add(jcbEnableAll);

        jcbFilterEnabled.setText(LABELS.LABEL_TABLE_HIDE_DISABLED);
        jcbFilterEnabled.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jcbFilterEnabledItemStateChanged(evt);
            }
        });
        jpCheckAll.add(jcbFilterEnabled);

        jpCenter.add(jpCheckAll, java.awt.BorderLayout.PAGE_START);
        jpCenter.add(jsPanel, java.awt.BorderLayout.CENTER);

        getContentPane().add(jpCenter, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void bSaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_bSaveActionPerformed
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
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton bSave;
    private javax.swing.JCheckBox jcbEnableAll;
    private javax.swing.JCheckBox jcbFilterEnabled;
    private javax.swing.JPanel jpCenter;
    private javax.swing.JPanel jpCheckAll;
    private javax.swing.JScrollPane jsPanel;
    private javax.swing.JToolBar toolBar;
    // End of variables declaration//GEN-END:variables

    class JPopUpMenuMicrograph extends JPopupMenu {

        private JMenuItem jmiShowCTF = new JMenuItem(LABELS.LABEL_SHOW_CTF);
        private JMenuItem jmiExtractColumn = new JMenuItem(LABELS.LABEL_EXTRACT_COLUMN);
        private JMenuItem jmiRecalculateCTF = new JMenuItem(LABELS.LABEL_RECALCULATE_CTF);

        public JPopUpMenuMicrograph() {
            super();

            add(jmiShowCTF);
            add(jmiExtractColumn);
            add(new JSeparator());
            add(jmiRecalculateCTF);

            jmiShowCTF.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    showCTFFile(table.getSelectedRow());
                }
            });

            jmiExtractColumn.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    extractColumn(table.getSelectedColumn());
                }
            });

            jmiRecalculateCTF.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    recalculateCTF(table.getSelectedRow());
                }
            });
        }

        private void refreshItems(int row, int col) {
            jmiShowCTF.setEnabled(tableModel.hasCtfData());
            jmiExtractColumn.setEnabled(table.getValueAt(row, col) instanceof TableImageItem);
        }

        private void showCTFFile(int row) {
            String CTFfile = tableModel.getCTFfile(row);

            frameCTF.loadFile(CTFfile);

            frameCTF.setLocationRelativeTo(this);
            frameCTF.setVisible(true);
        }

        private void extractColumn(int column) {
            LinkedList<String> filenames = new LinkedList<String>();

            for (int i = 0; i < table.getRowCount(); i++) {
                TableImageItem item = (TableImageItem) table.getValueAt(i, column);

                if (tableModel.isRowEnabled(i)) {
                    filenames.add(item.getFileName());
                }
            }

            ImagesWindowFactory.openTable(filenames.toArray(new String[filenames.size()]));
        }

        private void recalculateCTF(int row) {
            Object item = table.getValueAt(row, CTF_IMAGE_COLUMN);

            if (item instanceof TableImageItem) {
                showCTFImage((TableImageItem) item, tableModel.getCTFfile(row));
            }
        }
    }
}
