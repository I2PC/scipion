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

import browser.imageitems.TableImageItem;
import browser.table.renderers.DoubleRenderer;
import browser.table.renderers.FileNameRenderer;
import browser.table.renderers.ImageMicrographRenderer;
import browser.windows.ImagesWindowFactory;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameMicrographs extends JFrame {

    private final static String BUTTON_SAVE = "Save";
    private final static int CELL_WIDTH_MIN = 50;
    private final static int CELL_MARGIN = 2;
    private JTable table;
    private TableModelMicrographs tableModel;
    private XTableColumnModel columnModel = new XTableColumnModel();
    private JPopupMenu jPopupMenu = new JPopupMenu();
    private JMenuItem jmiShowCTF = new JMenuItem("Show CTF");
    private JFileChooser fc = new JFileChooser();
    FileNameRenderer fileNameRenderer = new FileNameRenderer();
    ImageMicrographRenderer imageRenderer = new ImageMicrographRenderer();
    DoubleRenderer doubleRenderer = new DoubleRenderer();
    private JFrameCTF frameCTF = new JFrameCTF();

    /** Creates new form JFrameMicrographs */
    public JFrameMicrographs(MetaData md) {
        super();

        setTitle(md.getFilename());

        initComponents();

        // Builds table.
        tableModel = new TableModelMicrographs(md);

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

        table.setAutoCreateRowSorter(true);
        table.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
        table.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                tableMouseClicked(evt);
            }
        });

        scrollPane.setViewportView(table);

        table.setColumnModel(columnModel);
        table.createDefaultColumnsFromModel();

        setRenderers();

        // places Defocus columns next to images.
        table.moveColumn(table.getColumnCount() - 2, TableModelMicrographs.DEFOCUS_U_COL);
        table.moveColumn(table.getColumnCount() - 1, TableModelMicrographs.DEFOCUS_V_COL);

        hideColumns();

        jPopupMenu.add(jmiShowCTF);
        jmiShowCTF.setEnabled(tableModel.hasCtfData());

        jmiShowCTF.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                showCTF(table.getSelectedRows()[0]);
            }
        });

        packColumns();
        packRows();

        pack();
    }

    private void hideColumns() {
        for (int i = 0; i < tableModel.columnsToHide.length; i++) {
            int index = tableModel.columnsToHide[i];
            TableColumn column = columnModel.getColumnByModelIndex(index);
            columnModel.setColumnVisible(column, false);
        }
    }

    private void showCTF(int row) {
        String CTFfile = tableModel.getCTFfile(row);

        frameCTF.setText(CTFfile);

        frameCTF.setLocationRelativeTo(this);
        frameCTF.setVisible(true);
    }

    private void setRenderers() {
        // Images
        setRenderer(table, TableModelMicrographs.imagesColumnIndex, imageRenderer);

        // Doubles
        setRenderer(table, TableModelMicrographs.doubleColumnIndex, doubleRenderer);

        // Filenames
        setRenderer(table, TableModelMicrographs.filenameColumnIndex, fileNameRenderer);
    }

    private static void setRenderer(JTable table, int indexes[], TableCellRenderer renderer) {
        for (int i = 0; i < indexes.length; i++) {
            if (table.getColumnCount() > indexes[i]) {
                table.getColumnModel().getColumn(indexes[i]).setCellRenderer(renderer);
            }
        }
    }

    private void packColumns() {
        for (int i = 0; i < table.getColumnCount(); i++) {
            packColumn(table, i, CELL_MARGIN);
        }
    }

    // Sets the preferred width of the visible column specified by vColIndex. The column
    // will be just wide enough to show the column head and the widest cell in the column.
    // margin pixels are added to the left and right
    // (resulting in an additional width of 2*margin pixels).
    private void packColumn(JTable table, int vColIndex, int margin) {
        DefaultTableColumnModel colModel = (DefaultTableColumnModel) table.getColumnModel();
        TableColumn col = colModel.getColumn(vColIndex);
        int width = 0;

        // Get width of column header
        TableCellRenderer renderer = col.getHeaderRenderer();
        if (renderer == null) {
            renderer = table.getTableHeader().getDefaultRenderer();
        }

//        Component comp = renderer.getTableCellRendererComponent(
//                table, col.getHeaderValue(), false, false, 0, 0);
        width = CELL_WIDTH_MIN;//comp.getPreferredSize().width;
        Component comp;

        // Get maximum width of column data
        for (int r = 0; r < table.getRowCount(); r++) {
            renderer = table.getCellRenderer(r, vColIndex);
            comp = renderer.getTableCellRendererComponent(
                    table, table.getValueAt(r, vColIndex), false, false, r, vColIndex);
            width = Math.max(width, comp.getPreferredSize().width);
        }

        // Add margin
        width += 2 * margin;

        // Set the width
        col.setPreferredWidth(width);
    }

    // The height of each row is set to the preferred height of the tallest cell in that row.
    public void packRows() {
        packRows(table, 0, table.getRowCount(), CELL_MARGIN);
    }

    // For each row >= start and < end, the height of a row is set to
    // the preferred height of the tallest cell in that row.
    public void packRows(JTable table, int start, int end, int margin) {
        for (int r = 0; r < table.getRowCount(); r++) {
            // Get the preferred height
            int h = getPreferredRowHeight(table, r, margin);

            // Now set the row height using the preferred height
            if (table.getRowHeight(r) != h) {
                table.setRowHeight(r, h);
            }
        }
    }

    private void tableMouseClicked(java.awt.event.MouseEvent evt) {
        int row = table.rowAtPoint(evt.getPoint());
        int col = table.columnAtPoint(evt.getPoint());

        if (SwingUtilities.isLeftMouseButton(evt)) {
            if (evt.getClickCount() > 1) {
                Object item = table.getValueAt(row, col);

                if (item instanceof TableImageItem) {
                    ImagesWindowFactory.openImage(((TableImageItem) item).getImagePlus());
                }
            } else if (col == 0) {  // Enabled / Disabled -> Update entire row.
                tableModel.fireTableRowsUpdated(row, row);
            }
        } else if (SwingUtilities.isRightMouseButton(evt)) {
            ListSelectionModel model = table.getSelectionModel();

            model.setSelectionInterval(row, row);

            jPopupMenu.show(evt.getComponent(), evt.getX(), evt.getY());
        }
    }

    // Returns the preferred height of a row.
    // The result is equal to the tallest cell in the row.
    public int getPreferredRowHeight(JTable table, int rowIndex, int margin) {
        int height = ImageMicrographRenderer.CELL_HEIGHT + 2 * margin;

        return height;
    }

    private void save(String fileName) {
        if (tableModel.save(table, fileName)) {
            JOptionPane.showMessageDialog(this, "File " + fileName + " sucessfully saved.", "File saved.", JOptionPane.INFORMATION_MESSAGE);
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

        scrollPane = new javax.swing.JScrollPane();
        toolBar = new javax.swing.JToolBar();
        bSave = new javax.swing.JButton();

        getContentPane().add(scrollPane, java.awt.BorderLayout.CENTER);

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

        getContentPane().add(toolBar, java.awt.BorderLayout.PAGE_START);

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
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton bSave;
    private javax.swing.JScrollPane scrollPane;
    private javax.swing.JToolBar toolBar;
    // End of variables declaration//GEN-END:variables
}
