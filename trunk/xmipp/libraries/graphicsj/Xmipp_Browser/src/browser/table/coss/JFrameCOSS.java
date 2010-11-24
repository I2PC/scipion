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

import browser.table.coss.renderers.ImageRenderer;
import browser.table.coss.renderers.InfoRenderer;
import java.awt.Component;
import java.awt.print.PrinterException;
import java.io.File;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameCOSS extends JFrame {

    private final static String BUTTON_SAVE = "Save";
    private final static String BUTTON_PRINT = "Print";
    private final static String MESSAGE_ERROR_PRINTER = "Printer Error";
    private final static int CELL_MARGIN = 2;
    private TableModelCOSS tableModel = new TableModelCOSS();
    private JFileChooser fc = new JFileChooser();
    ImageRenderer imageRenderer = new ImageRenderer();
    InfoRenderer infoRenderer = new InfoRenderer();
    ImageCellEditor imageCellEditor = new ImageCellEditor(this);

    /** Creates new form JFrameCOSS */
    public JFrameCOSS(String file) {
        initComponents();

        setTitle(file);

        tableModel.setRowSorter(table.getRowSorter());
        tableModel.setHeader(table.getTableHeader());

        String images[] = new String[]{
            "/home/juanjo/Desktop/SampleData/Imgs/img000001.xmp",
            "/home/juanjo/Desktop/SampleData/Imgs/img000002.xmp",
            "/home/juanjo/Desktop/SampleData/Imgs/img000003.xmp",
        };
        /*            "/home/juanjo/Desktop/Xmipp_Browser/src/resources/type_image.png",
        "/home/juanjo/Desktop/Xmipp_Browser/src/resources/type_folder.png",
        "/home/juanjo/Desktop/Xmipp_Browser/src/resources/capture.png"};*/

        for (int i = 0; i < 5; i++) {
            CossEntry entry = new CossEntry(i % 2 == 0,
                    images[i % images.length],
                    images[(i + 1) % images.length],
                    images[(i + 2) % images.length],
                    i + " info...");
            tableModel.addRow(new Object[]{new Boolean(entry.enabled), entry.img1, entry.img2, entry.img3, entry.info});
        }

        table.getColumnModel().getColumn(1).setCellRenderer(imageRenderer);
        table.getColumnModel().getColumn(2).setCellRenderer(imageRenderer);
        table.getColumnModel().getColumn(3).setCellRenderer(imageRenderer);
        table.getColumnModel().getColumn(4).setCellRenderer(infoRenderer);

        table.getColumnModel().getColumn(1).setCellEditor(imageCellEditor);
        table.getColumnModel().getColumn(2).setCellEditor(imageCellEditor);
        table.getColumnModel().getColumn(3).setCellEditor(imageCellEditor);

        packRows();
    }

    // The height of each row is set to the preferred height of the tallest cell in that row.
    public void packRows() {
        packRows(table, 0, table.getRowCount(), CELL_MARGIN);
    }

    // Returns the preferred height of a row.
    // The result is equal to the tallest cell in the row.
    public int getPreferredRowHeight(JTable table, int rowIndex, int margin) {
        // Get the current default height for all rows
        int height = table.getRowHeight();

        // Determine highest cell in the row
        for (int c = 0; c < table.getColumnCount(); c++) {
            TableCellRenderer renderer = table.getCellRenderer(rowIndex, c);
            Component comp = table.prepareRenderer(renderer, rowIndex, c);
            int h = comp.getPreferredSize().height + 2 * margin;
            height = Math.max(height, h);
        }
        return height;
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

    private void save(File file) {
        tableModel.save(file);
    }

    private void removeRows() {
        int rows[] = table.getSelectedRows();

        for (int i = rows.length; i > 0; i--) {
            tableModel.removeRow(rows[i - 1]);
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
        bAddRows = new javax.swing.JButton();
        bRemoveRows = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        table.setAutoCreateRowSorter(true);
        table.setModel(tableModel);
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

        bAddRows.setText("+Row");
        bAddRows.setFocusable(false);
        bAddRows.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        bAddRows.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        bAddRows.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                bAddRowsActionPerformed(evt);
            }
        });
        toolBar.add(bAddRows);

        bRemoveRows.setText("-Row(s)");
        bRemoveRows.setFocusable(false);
        bRemoveRows.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        bRemoveRows.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        bRemoveRows.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                bRemoveRowsActionPerformed(evt);
            }
        });
        toolBar.add(bRemoveRows);

        getContentPane().add(toolBar, java.awt.BorderLayout.PAGE_START);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void bSaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_bSaveActionPerformed
        if (fc.showSaveDialog(this) != JFileChooser.CANCEL_OPTION) {
            save(fc.getSelectedFile());
        }
    }//GEN-LAST:event_bSaveActionPerformed

    private void bPrintActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_bPrintActionPerformed
        printTable();
    }//GEN-LAST:event_bPrintActionPerformed

    private void bAddRowsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_bAddRowsActionPerformed
        tableModel.appendEmptyRow();
    }//GEN-LAST:event_bAddRowsActionPerformed

    private void bRemoveRowsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_bRemoveRowsActionPerformed
        removeRows();
    }//GEN-LAST:event_bRemoveRowsActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton bAddRows;
    private javax.swing.JButton bPrint;
    private javax.swing.JButton bRemoveRows;
    private javax.swing.JButton bSave;
    private javax.swing.JScrollPane scrollPane;
    private javax.swing.JTable table;
    private javax.swing.JToolBar toolBar;
    // End of variables declaration//GEN-END:variables
}
