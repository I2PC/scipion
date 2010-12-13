/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * JFrameCOSS.java
 *
 * Created on Oct 25, 2010, 10:50:01 AM
 */
package browser.table.coss;

import browser.table.coss.renderers.FileNameRenderer;
import browser.table.coss.renderers.ImageRenderer;
import browser.table.coss.renderers.InfoRenderer;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.print.PrinterException;
import java.io.File;
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
import javax.swing.table.TableColumnModel;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameCOSS extends JFrame {

    private final static String BUTTON_SAVE = "Save";
    private final static String BUTTON_PRINT = "Print";
    private final static String MESSAGE_ERROR_PRINTER = "Printer Error";
    private final static int CELL_MARGIN = 2;
    private final static int CELL_WIDTH_MIN = 50;
    private TableModelCOSS tableModel;
    private XTableColumnModel columnModel = new XTableColumnModel();
    private JPopupMenu jPopupMenu = new JPopupMenu();
    private JMenuItem jmiShowCTF = new JMenuItem("Show CTF");
    private JFileChooser fc = new JFileChooser();
    ImageRenderer imageRenderer = new ImageRenderer();
    FileNameRenderer fileNameRenderer = new FileNameRenderer();
    InfoRenderer infoRenderer = new InfoRenderer();
    JFrameCTF frameCTF = new JFrameCTF();
    _ImageCellEditor imageCellEditor = new _ImageCellEditor(this);
    String baseDir;

    static {
        System.loadLibrary("XmippDataJava");
    }

    /** Creates new form JFrameCOSS */
    public JFrameCOSS(String fileName) {
        super(fileName);

        baseDir = new File(fileName).getParent();  // Stores base dir.

        // Builds table.
        tableModel = new TableModelCOSS(fileName);

        initComponents();

        table.setColumnModel(columnModel);
        table.createDefaultColumnsFromModel();

        setInvisibleColumns();
        setRenderers();

        jPopupMenu.add(jmiShowCTF);

        jmiShowCTF.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                showCTF(table.getSelectedRows()[0]);
            }
        });

        packColumns();
        packRows();
    }

    private void setInvisibleColumns() {
        for (int i = 0; i < TableModelCOSS.invisibleColumns.length; i++) {
            int index = TableModelCOSS.invisibleColumns[i];
            TableColumn column = columnModel.getColumnByModelIndex(index);
            columnModel.setColumnVisible(column, false);
        }
    }

    private void showCTF(int row) {
        String CTFfile = tableModel.getCTFfile(table, row);

        frameCTF.setText(baseDir, CTFfile);

        frameCTF.setLocationRelativeTo(this);
        frameCTF.setVisible(true);
    }

    private void setRenderers() {
        for (int i = 0; i < TableModelCOSS.imagesColumnIndex.length; i++) {
            int j = TableModelCOSS.imagesColumnIndex[i];
            table.getColumnModel().getColumn(j).setCellRenderer(imageRenderer);
        }

        table.getColumnModel().getColumn(1).setCellRenderer(fileNameRenderer);
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

    // Returns the preferred height of a row.
    // The result is equal to the tallest cell in the row.
    public int getPreferredRowHeight(JTable table, int rowIndex, int margin) {
        // Get the current default height for all rows
//        int height = table.getRowHeight();
//
//        // Determine highest cell in the row
//        for (int c = 0; c < table.getColumnCount(); c++) {
//            System.out.println(" ####### Column: " + c);
//            TableCellRenderer renderer = table.getCellRenderer(rowIndex, c);
//            Component comp = table.prepareRenderer(renderer, rowIndex, c);
//            int h = comp.getPreferredSize().height + 2 * margin;
//            height = Math.max(height, h);
//        }

        int height = TableImageItemCOSS.CELL_HEIGHT + 2 * margin;

        return height;
    }

    private void save(String fileName) {
        if (TableModelCOSS.save(table, fileName)) {
            JOptionPane.showMessageDialog(this, "File " + fileName + " sucessfully saved.", "File saved.", JOptionPane.INFORMATION_MESSAGE);
        }
    }

    private void printTable() {
        try {
            table.print();
        } catch (PrinterException pex) {
            pex.printStackTrace();
            JOptionPane.showMessageDialog(this, pex.getMessage(), MESSAGE_ERROR_PRINTER, JOptionPane.ERROR_MESSAGE);
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
        table = new javax.swing.JTable();
        toolBar = new javax.swing.JToolBar();
        bSave = new javax.swing.JButton();
        bPrint = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        table.setAutoCreateRowSorter(true);
        table.setModel(tableModel);
        table.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
        table.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                tableMouseClicked(evt);
            }
        });
        scrollPane.setViewportView(table);

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

        bPrint.setText(BUTTON_PRINT);
        bPrint.setFocusable(false);
        bPrint.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        bPrint.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        bPrint.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                bPrintActionPerformed(evt);
            }
        });
        toolBar.add(bPrint);

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

    private void bPrintActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_bPrintActionPerformed
        printTable();
    }//GEN-LAST:event_bPrintActionPerformed

    private void tableMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_tableMouseClicked
        if (SwingUtilities.isLeftMouseButton(evt)) {
            if (evt.getClickCount() > 1) {
                int rows[] = table.getSelectedRows();
                int cols[] = table.getSelectedColumns();

                Object item = table.getValueAt(
                        rows[0], cols[0]);

                if (item instanceof TableImageItemCOSS) {
                    ((TableImageItemCOSS) item).getFullZiseImagePlus().show();
                }
            }
        } else if (SwingUtilities.isRightMouseButton(evt)) {
            int rowNumber = table.rowAtPoint(evt.getPoint());

            ListSelectionModel model = table.getSelectionModel();

            model.setSelectionInterval(rowNumber, rowNumber);

            jPopupMenu.show(evt.getComponent(), evt.getX(), evt.getY());
        }
    }//GEN-LAST:event_tableMouseClicked
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton bPrint;
    private javax.swing.JButton bSave;
    private javax.swing.JScrollPane scrollPane;
    private javax.swing.JTable table;
    private javax.swing.JToolBar toolBar;
    // End of variables declaration//GEN-END:variables

    class ToolTipHeader extends JTableHeader {

        String[] toolTips;

        public ToolTipHeader(TableColumnModel model) {
            super(model);
        }

        @Override
        public String getToolTipText(MouseEvent e) {
            int col = columnAtPoint(e.getPoint());
            int modelCol = getTable().convertColumnIndexToModel(col);
            String retStr;
            try {
                retStr = toolTips[modelCol];
            } catch (NullPointerException ex) {
                retStr = "";
            } catch (ArrayIndexOutOfBoundsException ex) {
                retStr = "";
            }
            if (retStr.length() < 1) {
                retStr = super.getToolTipText(e);
            }
            return retStr;
        }

        public void setToolTipStrings(String[] toolTips) {
            this.toolTips = toolTips;
        }
    }
}
