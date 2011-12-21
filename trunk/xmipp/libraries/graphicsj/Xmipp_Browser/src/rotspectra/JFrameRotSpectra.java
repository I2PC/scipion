/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * JFrameRotSpectra.java
 *
 * Created on Sep 23, 2011, 4:31:29 PM
 */
package rotspectra;

import browser.DEBUG;
import browser.LABELS;
import browser.gallery.JFrameGallery;
import browser.gallery.models.GalleryRowHeaderModel;
import browser.gallery.renderers.RowHeaderRenderer;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Enumeration;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.LookAndFeel;
import javax.swing.SwingUtilities;
import javax.swing.border.Border;
import javax.swing.table.TableColumn;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameRotSpectra extends javax.swing.JFrame {

    final static int DEFAULT_WIDTH = 128;
    final static int DEFAULT_HEIGHT = 128;
    private RotSpectraTableModel tableModel;
    //private ImagesTableColumnModel columnModel;
    private GalleryRowHeaderModel rowHeaderModel;
    private int previousSelectedRow, previousSelectedCol;
    private JList rowHeader;
    private RotSpectraRenderer renderer = new RotSpectraRenderer(DEFAULT_WIDTH, DEFAULT_HEIGHT);
    private JPopUpMenuTable jpopUpMenuTable;
//
//    public static void main(String args[]) {
//        try {
//            String dir = "/data2/MicrographPreprocessing/2D/RotSpectra/run_001/";
//            String filenameClasses = dir + "results_classes.xmd";
//            String filenameVectors = dir + "results_vectors.xmd";
//            String vectorsData = dir + "results_vectors.xmd.raw";
//
//            JFrameRotSpectra frame = new JFrameRotSpectra(filenameVectors, filenameClasses, vectorsData);
//            ImagesWindowFactory.setConvenientSize(frame);
//            frame.setLocationRelativeTo(null);
//            frame.setVisible(true);
//        } catch (Exception ex) {
//            DEBUG.printException(ex);
//        }
//    }

    public JFrameRotSpectra(String fnVectors, String fnClasses, String fnData) {
        super();

        try {
            tableModel = new RotSpectraTableModel(fnClasses, fnVectors, fnData);

            initComponents();

            setTitle(fnVectors + ": " + tableModel.getSize() + " vectors.");

            table.setDefaultRenderer(RotSpectraVector.class, renderer);

            setRowHeader();

            packColumns();
            packRows();

            jpopUpMenuTable = new JPopUpMenuTable();
        } catch (Exception e) {
            DEBUG.printException(e);
            IJ.error(e.getMessage());
        }
    }

    private void setRowHeader() {
        rowHeaderModel = new GalleryRowHeaderModel(table, 1);

        rowHeader = new JList();
        rowHeader.setModel(rowHeaderModel);

        LookAndFeel.installColorsAndFont(rowHeader, "TableHeader.background",
                "TableHeader.foreground", "TableHeader.font");

        rowHeader.setCellRenderer(new RowHeaderRenderer());

        jspTable.setRowHeaderView(rowHeader);
//        rowHeader.repaint();
    }

    private void packRows() {
        //  Row header.
        rowHeader.setFixedCellHeight(renderer.getHeight());

        for (int row = 0; row < table.getRowCount(); row++) {
            table.setRowHeight(row, renderer.getHeight());
        }
    }

    private void packColumns() {
        Enumeration<TableColumn> columns = table.getColumnModel().getColumns();
        while (columns.hasMoreElements()) {
            columns.nextElement().setPreferredWidth(renderer.getWidth());
        }
    }

    private void setShowLabels(boolean show) {
        tableModel.setShowLabels(show);

        updateTable();
    }

    private synchronized void updateTable() {
        if (tableModel.getSize() > 0) {
            Dimension dimension = getCellSize();
            renderer.setPreferredSize(dimension);

            // Adjusts rows size.
            table.setRowHeight(dimension.height > 0 ? dimension.height : 1);
            rowHeader.setFixedCellHeight(table.getRowHeight());
//            rowHeader.setFixedCellWidth(rowHeader.getModel().getSize() * 2);
//            columnModel.setWidth(dimension.width > 0 ? dimension.width : 1);

            rowHeader.repaint();
            table.revalidate();
        }
    }

    private Dimension getCellSize() {
        int font_height = 0;
        if (tableModel.isShowingLabels()) {
            font_height = renderer.getFontMetrics(renderer.getFont()).getHeight();
            font_height += renderer.getIconTextGap();  // Adds the extra gap.
            font_height -= table.getIntercellSpacing().height;   // Removes spacing.
        }

        int borderHeight = 0, borderWidth = 0;
        Border border = renderer.getBorder();

        if (border != null) {
            Insets insets = renderer.getBorder().getBorderInsets(renderer);
            borderWidth = insets.left + insets.right;
            borderHeight = insets.bottom + insets.top;
        }

        return new Dimension(
                DEFAULT_WIDTH + 2 * borderWidth,
                DEFAULT_HEIGHT + font_height + 2 * borderHeight);
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jpControls = new javax.swing.JPanel();
        jcbShowLabels = new javax.swing.JCheckBox();
        jspTable = new javax.swing.JScrollPane();
        table = new javax.swing.JTable();

        jpControls.setLayout(new javax.swing.BoxLayout(jpControls, javax.swing.BoxLayout.LINE_AXIS));

        jcbShowLabels.setText(LABELS.LABEL_ROTSPECTRA_SHOW_NIMAGES);
        jcbShowLabels.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jcbShowLabelsActionPerformed(evt);
            }
        });
        jpControls.add(jcbShowLabels);

        getContentPane().add(jpControls, java.awt.BorderLayout.NORTH);

        table.setModel(tableModel);
        table.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
        table.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                tableMouseClicked(evt);
            }
        });
        jspTable.setViewportView(table);

        getContentPane().add(jspTable, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void tableMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_tableMouseClicked
        int view_row = table.rowAtPoint(evt.getPoint());
        int view_col = table.columnAtPoint(evt.getPoint());

        if (evt.getButton() == MouseEvent.BUTTON1) {  // Left click.
            if (evt.getClickCount() > 1) {
                RotSpectraVector rsVector = (RotSpectraVector) table.getValueAt(view_row, view_col);

                JFrame frameChart = rsVector.getChart();
                frameChart.setVisible(true);
            } else {
                // Ctrl adds items to selection, otherwise previous ones are removed.
                if (!evt.isControlDown()) {
                    tableModel.clearSelection();
                }

                if (evt.isShiftDown()) {
                    tableModel.selectRange(
                            previousSelectedRow, previousSelectedCol,
                            view_row, view_col);
                } else {
                    tableModel.setSelected(view_row, view_col, true);
                }

                table.repaint();
            }

            if (!evt.isShiftDown()) {
                previousSelectedRow = view_row;
                previousSelectedCol = view_col;
            }
        } else if (evt.getButton() == MouseEvent.BUTTON3) {  // Right click.
            if (tableModel.getSelectedItems().size() < 2) {
                tableModel.clearSelection();
                tableModel.setSelected(view_row, view_col, true);

                table.setRowSelectionInterval(view_row, view_row);
                table.setColumnSelectionInterval(view_col, view_col);

                table.repaint();
            }

            final MouseEvent me = evt;
            SwingUtilities.invokeLater(new Runnable() {

                public void run() {
                    jpopUpMenuTable.show(me.getComponent(), me.getPoint());
                }
            });
        }
    }//GEN-LAST:event_tableMouseClicked

    private void jcbShowLabelsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jcbShowLabelsActionPerformed
        setShowLabels(jcbShowLabels.isSelected());
    }//GEN-LAST:event_jcbShowLabelsActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JCheckBox jcbShowLabels;
    private javax.swing.JPanel jpControls;
    private javax.swing.JScrollPane jspTable;
    private javax.swing.JTable table;
    // End of variables declaration//GEN-END:variables

    class JPopUpMenuTable extends JPopupMenu {

        Point location;
        JMenuItem jmiOpenAsGallery = new JMenuItem("Open as Gallery");

        public JPopUpMenuTable() {
            add(jmiOpenAsGallery);

            jmiOpenAsGallery.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
//                    int row = table.rowAtPoint(location);
//                    int col = table.columnAtPoint(location);

                    // Gets all images contained by selected vectors.
                    ArrayList<String> filenames = new ArrayList<String>();
                    for (int i = 0; i < tableModel.getSelectedItems().size(); i++) {
                        RotSpectraVector vector = tableModel.getSelectedItems().get(i);
                        filenames.addAll(vector.getImagesFilenames());
                    }

                    String array[] = filenames.toArray(new String[filenames.size()]);

                    JFrameGallery frame = ImagesWindowFactory.openFilesAsGallery(array, true);
                    frame.setTitle(tableModel.getSelectedItems().get(0).block
                            + "...: " + filenames.size() + " images.");
                }
            });
        }

        public void show(Component cmpnt, Point location) {
            this.location = location;

            int row = table.rowAtPoint(location);
            int col = table.columnAtPoint(location);
            RotSpectraVector rsv = (RotSpectraVector) table.getValueAt(row, col);

            jmiOpenAsGallery.setEnabled(rsv.getNImages() > 0);

            show(cmpnt, location.x, location.y);
        }
    }
}
