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

import browser.LABELS;
import browser.imageitems.TableImageItem;
import browser.table.ImagesRowHeaderModel;
import browser.table.micrographs.filters.EnableFilter;
import browser.table.micrographs.renderers.MicrographDoubleRenderer;
import browser.table.micrographs.renderers.MicrographFileNameRenderer;
import browser.table.micrographs.renderers.MicrographImageRenderer;
import browser.table.micrographs.singlecolumntable.JFrameExtractColumn;
import browser.table.renderers.RowHeaderRenderer;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import ij.gui.GenericDialog;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.Vector;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
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
public class JFrameMicrographs extends JFrame {

    private final static String BUTTON_SAVE = "Save";
    private final static int CELL_WIDTH_MIN = 50;
    private JTable table;
    private MicrographsTableModel tableModel;
    private ImagesRowHeaderModel rowHeaderModel;
    private XTableColumnModel columnModel = new XTableColumnModel();
    private JPopupMenu jPopupMenu = new JPopupMenu();
    private JMenuItem jmiShowCTF = new JMenuItem("Show CTF");
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
        table.moveColumn(table.getColumnCount() - 2, MicrographsTableModel.DEFOCUS_U_COL);
        table.moveColumn(table.getColumnCount() - 1, MicrographsTableModel.DEFOCUS_V_COL);

        hideColumns();
        setRowHeader();

        jPopupMenu.add(jmiShowCTF);
        jmiShowCTF.setEnabled(tableModel.hasCtfData());

        jmiShowCTF.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                showCTF(table.getSelectedRows()[0]);
            }
        });

        updateTable();

        pack();
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

    private void showCTF(int row) {
        String CTFfile = tableModel.getCTFfile(row);

        frameCTF.loadFile(CTFfile);

        frameCTF.setLocationRelativeTo(this);
        frameCTF.setVisible(true);
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
        // Images
        setRenderer(table, MicrographsTableModel.imagesColumnIndex, imageRenderer);

        // Doubles
        setRenderer(table, MicrographsTableModel.doubleColumnIndex, doubleRenderer);

        // Filenames
        setRenderer(table, MicrographsTableModel.filenameColumnIndex, fileNameRenderer);
    }

    private static void setRenderer(JTable table, int indexes[], TableCellRenderer renderer) {
        for (int i = 0; i < indexes.length; i++) {
            if (table.getColumnCount() > indexes[i]) {
                table.getColumnModel().getColumn(indexes[i]).setCellRenderer(renderer);
            }
        }
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
    private void packColumn(JTable table, int vColIndex) {
        DefaultTableColumnModel colModel = (DefaultTableColumnModel) table.getColumnModel();
        TableColumn col = colModel.getColumn(vColIndex);
        int width = 0;

        // Get width of column header
        TableCellRenderer renderer = col.getHeaderRenderer();
        if (renderer == null) {
            renderer = table.getTableHeader().getDefaultRenderer();
        }

        width = CELL_WIDTH_MIN;

        // Get maximum width of column data
        for (int r = 0; r < table.getRowCount(); r++) {
            renderer = table.getCellRenderer(r, vColIndex);
            Component comp = renderer.getTableCellRendererComponent(
                    table, table.getValueAt(r, vColIndex), false, false, r, vColIndex);
            width = Math.max(width, comp.getPreferredSize().width);
        }

        // Set the width
        col.setPreferredWidth(width);
    }

    // The height of each row is set to the preferred height of the tallest cell in that row.
    private void packRows() {
        packRows(table, 0, table.getRowCount());
    }

    // For each row >= start and < end, the height of a row is set to
    // the preferred height of the tallest cell in that row.
    private void packRows(JTable table, int start, int end) {
        for (int r = 0; r < table.getRowCount(); r++) {

            // Now set the row height using the preferred height
            if (table.getRowHeight(r) != MicrographImageRenderer.CELL_HEIGHT) {
                table.setRowHeight(r, MicrographImageRenderer.CELL_HEIGHT);
            }
        }
    }

    private void tableMouseClicked(java.awt.event.MouseEvent evt) {
        int row = table.convertRowIndexToModel(table.rowAtPoint(evt.getPoint()));
        int col = table.convertColumnIndexToModel(table.columnAtPoint(evt.getPoint()));

        if (SwingUtilities.isLeftMouseButton(evt)) {
            if (evt.getClickCount() > 1) {
                Object item = table.getValueAt(row, col);

                if (item instanceof TableImageItem) {
                    ImagesWindowFactory.openImage(((TableImageItem) item).getImagePlus());
                }
            } else {    // Single click
                int min = table.getSelectionModel().getMinSelectionIndex();
                int max = table.getSelectionModel().getMaxSelectionIndex();

                if (evt.isShiftDown()) {
                    // Last selection value will be the value for all the selected items.
                    Boolean value = (Boolean) table.getValueAt(row, 0);

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
            table.getSelectionModel().setSelectionInterval(row, row);

            jPopupMenu.show(evt.getComponent(), evt.getX(), evt.getY());
        }
    }

    public static void main(String args[]) {
        JFrameMicrographs frame = new JFrameMicrographs("/media/PENDRIVE/Ad5GLflagIIIa/Preprocessing/all_micrographs_en_test.sel");
        frame.setVisible(true);
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

    private void extractColumn() {
        int index = -1;
        final GenericDialog dialog = new GenericDialog("Select column to extract: ");
        final String[] columns = new String[table.getColumnCount()];

        // Builds choice.
        for (int i = 0; i < columns.length; i++) {
            int column = table.convertColumnIndexToModel(i);
            columns[i] = tableModel.getColumnName(column);
        }
        dialog.addChoice("column:", columns, columns[0]);
        dialog.showDialog();

        if (!dialog.wasCanceled()) {
            index = dialog.getNextChoiceIndex();
        }

        // If not cancelled...
        if (index >= 0) {
            JFrameExtractColumn frame = new JFrameExtractColumn(tableModel,
                    table.convertColumnIndexToModel(index));
            frame.pack();
            frame.setLocationRelativeTo(this);
            frame.setVisible(true); // ...shows it.
        }

        /*
        if (item.isStack()) {

        indexes[0] = "All";
        for (int i = ImageDouble.FIRST_IMAGE; i <= item.getNImages(); i++) {
        indexes[i] = String.valueOf(i);
        }

        dialog.addChoice("image:", indexes, indexes[0]);
        dialog.showDialog();
        if (!dialog.wasCanceled()) {
        image = dialog.getNextChoiceIndex();
        }
        } else {
        image = 0;
        }

        return image;
         */
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
        jbPrint = new javax.swing.JButton();
        jbExtractColumn = new javax.swing.JButton();
        jpCenter = new javax.swing.JPanel();
        jpCheckAll = new javax.swing.JPanel();
        jcbEnableAll = new javax.swing.JCheckBox();
        jcbFilterEnabled = new javax.swing.JCheckBox();
        jsPanel = new javax.swing.JScrollPane();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        toolBar.setRollover(true);

        bSave.setText(BUTTON_SAVE);
        bSave.setFocusable(false);
        bSave.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        bSave.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        bSave.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                bSaveActionPerformed(evt);
            }
        });
        toolBar.add(bSave);

        jbPrint.setText("[[Print]]");
        jbPrint.setFocusable(false);
        jbPrint.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbPrint.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbPrint.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbPrintActionPerformed(evt);
            }
        });
        toolBar.add(jbPrint);

        jbExtractColumn.setText("[[ExtractColumn]]");
        jbExtractColumn.setFocusable(false);
        jbExtractColumn.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbExtractColumn.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbExtractColumn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbExtractColumnActionPerformed(evt);
            }
        });
        toolBar.add(jbExtractColumn);

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

    private void jbPrintActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbPrintActionPerformed
        tableModel.print();
    }//GEN-LAST:event_jbPrintActionPerformed

    private void jbExtractColumnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbExtractColumnActionPerformed
        extractColumn();
    }//GEN-LAST:event_jbExtractColumnActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton bSave;
    private javax.swing.JButton jbExtractColumn;
    private javax.swing.JButton jbPrint;
    private javax.swing.JCheckBox jcbEnableAll;
    private javax.swing.JCheckBox jcbFilterEnabled;
    private javax.swing.JPanel jpCenter;
    private javax.swing.JPanel jpCheckAll;
    private javax.swing.JScrollPane jsPanel;
    private javax.swing.JToolBar toolBar;
    // End of variables declaration//GEN-END:variables
}
