/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * JFrameImagesTable.java
 *
 * Created on 22-abr-2010, 12:35:07
 */
package browser.table;

import browser.table.models.ImagesRowHeaderModel;
import browser.table.models.MDTableModel;
import browser.table.models.ImagesTableColumnModel;
import browser.imageitems.tableitems.AbstractTableImageItem;
import browser.table.renderers.ImageRenderer;
import browser.table.renderers.RowHeaderRenderer;
import browser.windows.ImagesWindowFactory;
import browser.DEBUG;
import browser.LABELS;
import browser.SpringUtilities;
import browser.imageitems.ImageConverter;
import browser.table.models.AbstractXmippTableModel;
import browser.table.models.VolumeTableModel;
import ij.IJ;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.ItemEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.util.Timer;
import java.util.TimerTask;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.KeyStroke;
import javax.swing.LookAndFeel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SpringLayout;
import javax.swing.SwingUtilities;
import javax.swing.border.Border;
import xmipp.Filename;
import xmipp.ImageDouble;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameImagesTable extends JFrame {//implements TableModelListener {

    private final static int DELAY_TO_UPDATE = 500;
    private AbstractXmippTableModel tableModel;
    private ImagesTableColumnModel columnModel;
    private ImagesRowHeaderModel rowHeaderModel;
    private int previousSelectedRow, previousSelectedCol;
    private JList rowHeader;
    private ImageRenderer renderer = new ImageRenderer();
    private Timer updateTimer = new Timer(true);    // Timer for zoom.
    private TableUpdater tableUpdaterTask;    // Associated task for zoom timer.
    private boolean isUpdating;
    private boolean autoAdjustColumns = false;
    private JPopUpMenuTable jpopUpMenuTable;
    private JMenuBarTable jMenuBarTable;
    JFileChooser fc = new JFileChooser();

    public JFrameImagesTable(String filename) {
        super();

        try {
            tableModel = Filename.isVolume(filename)
                    ? new VolumeTableModel(filename) : new MDTableModel(filename);

            postInit();

            // Set comboBox items.
            String labels[] = tableModel.getLabels();
            for (int i = 0; i < labels.length; i++) {
                jcbMDLabels.addItem(labels[i]);
            }
        } catch (Exception e) {
            DEBUG.printException(e);
            IJ.error(e.getMessage());
        }
    }

    public JFrameImagesTable(String filenames[]) {
        this(filenames, null);
    }

    public JFrameImagesTable(String filenames[], boolean enabled[]) {
        super();

        tableModel = new MDTableModel(filenames, enabled);
        postInit();
    }

    private void postInit() {
        jMenuBarTable = new JMenuBarTable();
        setJMenuBar(jMenuBarTable);

        columnModel = new ImagesTableColumnModel();

        initComponents();

        jtbUseGeometry.setEnabled(tableModel.isMetaData());    // Not aplicable for volumes.

        table.setColumnModel(columnModel);
        table.setDefaultRenderer(AbstractTableImageItem.class, renderer);

        setTitle(tableModel.getTitle());

        jpopUpMenuTable = new JPopUpMenuTable();

        // Sets limits for spinners.
        jsRows.setModel(new SpinnerNumberModel(1, 1, tableModel.getSize(), 1));
        jsColumns.setModel(new SpinnerNumberModel(1, 1, tableModel.getSize(), 1));
        jsGoToImage.setModel(new SpinnerNumberModel(1, 1, tableModel.getSize(), 1));

        setRowHeader();

        ((JSpinner.NumberEditor) jsZoom.getEditor()).getTextField().setColumns(8);
        ((JSpinner.NumberEditor) jsGoToImage.getEditor()).getTextField().setColumns(8);
        ((JSpinner.NumberEditor) jsRows.getEditor()).getTextField().setColumns(8);
        ((JSpinner.NumberEditor) jsColumns.getEditor()).getTextField().setColumns(8);

        jpDisplay.setLayout(new SpringLayout());
        SpringUtilities.makeCompactGrid(jpDisplay, 1, 7, 0, 0, 3, 3);
        jpStructure.setLayout(new SpringLayout());
        SpringUtilities.makeCompactGrid(jpStructure, 1, 5, 0, 0, 3, 3);

        // Stacks will be "auto-normalized".
        setNormalized(tableModel.isVolume());

        // Geometry info is used if present, otherwise button is disabled.
        boolean containsGeometry = tableModel.containsGeometryInfo();
        setUseGeometry(containsGeometry);
        jtbUseGeometry.setSelected(containsGeometry);
        jtbUseGeometry.setEnabled(containsGeometry);
    }

    private void setRowHeader() {
        rowHeaderModel = new ImagesRowHeaderModel(table, 1);

        rowHeader = new JList();
        rowHeader.setModel(rowHeaderModel);

        LookAndFeel.installColorsAndFont(rowHeader, "TableHeader.background",
                "TableHeader.foreground", "TableHeader.font");

        rowHeader.setCellRenderer(new RowHeaderRenderer());

        jsPanel.setRowHeaderView(rowHeader);
        rowHeader.repaint();
    }

    private void updateTable() {
//        if (table.isShowing()) {
        isUpdating = true;
        DEBUG.printMessage(" *** Updating table: " + System.currentTimeMillis());

        tableModel.updateSort();

        if (tableModel.getSize() > 0) {
            Dimension dimension = getCellSize();
            renderer.setPreferredSize(dimension);

            // Adjusts rows size.
            table.setRowHeight(dimension.height > 0 ? dimension.height : 1);
            rowHeader.setFixedCellHeight(table.getRowHeight());
//            rowHeader.setFixedCellWidth(rowHeader.getModel().getSize() * 2);
            columnModel.setWidth(dimension.width > 0 ? dimension.width : 1);

            // If auto adjust columns is enabled, refresh!
            if (autoAdjustColumns) {
                autoAdjustColumns();
            }

            rowHeader.repaint();
            table.revalidate();
        }

        isUpdating = false;
//        }
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
                tableModel.getCellWidth() + 2 * borderWidth,
                tableModel.getCellHeight() + font_height + 2 * borderHeight);
    }

    private void startUpdater() {
        if (tableUpdaterTask != null) {
            tableUpdaterTask.cancel();
        }

        tableUpdaterTask = new TableUpdater();
        updateTimer.schedule(tableUpdaterTask, DELAY_TO_UPDATE);
    }

    private void setShowLabels(boolean show) {
        tableModel.setShowLabels(show);

//        startUpdater();
        updateTable();
    }

    private void updateLabelToShow(int index) {
        tableModel.setSelectedLabel(index);

        updateTable();
    }

//    public void setDimensions(int rows, int columns) {
//        if (rows > 0) {
//            setRows(rows);            
//        }
//
//        if (columns > 0) {
//            setColumns(columns);
//        }
//    }
    private void autoAdjustColumns() {
        //System.out.println("rowHeader.getWidth(): " + rowHeader.getWidth());
        tableModel.autoAdjustColumns(
                //jsPanel.getVisibleRect().width - rowHeader.getWidth(),
                //jsPanel.getViewportBorderBounds().width - rowHeader.getWidth(),
                jsPanel.getViewport().getWidth() - rowHeader.getWidth(),
                table.getIntercellSpacing().width);

        jsRows.setValue(tableModel.getRowCount());
        jsColumns.setValue(tableModel.getColumnCount());
    }

    private void setInitialValues() {
        if (tableModel.getSize() > 0) {
            double scale = tableModel.getInitialZoomScale(
                    jsPanel.getVisibleRect().width,
                    table.getIntercellSpacing().width);
            setZoom((int) (scale * 100));
        } else {
            jsZoom.setEnabled(false);
            jcbAutoAdjustColumns.setEnabled(false);
            jcbShowLabels.setEnabled(false);
            jsRows.setEnabled(false);
            jsColumns.setEnabled(false);
            jsGoToImage.setEnabled(false);
        }

        updateTable();
    }

    private void setZoom(int zoom) {
        System.out.println(" *** Setting scale to: " + zoom / 100.0);
        isUpdating = true;
        tableModel.setZoomScale(zoom / 100.0);

        jsZoom.setValue(zoom);

        startUpdater();
        isUpdating = false;
    }

    private void setRows(int rows) {
        isUpdating = true;

        setRowsValue(rows);

        startUpdater();
        isUpdating = false;
    }

    private void setColumns(int columns) {
        isUpdating = true;

        setColumnsValue(columns);

        startUpdater();
        isUpdating = false;
    }

    public void setDimensions(int rows, int columns) {
        isUpdating = true;

        setAutoAdjustColumns(false);

        setRowsValue(rows);
        setColumnsValue(columns);

        startUpdater();
        isUpdating = false;
    }

    void setRowsValue(int rows) {
        if (rows > 0) {
            tableModel.setRows(rows);
            jsColumns.setValue(tableModel.getColumnCount());
        }
    }

    void setColumnsValue(int columns) {
        if (columns > 0) {
            tableModel.setColumns(columns);
            jsRows.setValue(tableModel.getRowCount());
        }
    }

    public void enableAutoadjustColumns(boolean enable) {
        jcbAutoAdjustColumns.setSelected(enable);
    }

    public void setAutoAdjustColumns(boolean autoAdjustColumns) {
//        if (!isUpdating) {
        isUpdating = true;

        this.autoAdjustColumns = autoAdjustColumns;

        jcbAutoAdjustColumns.setSelected(autoAdjustColumns);
        jsColumns.setEnabled(!autoAdjustColumns);
        jsRows.setEnabled(!autoAdjustColumns);

        if (autoAdjustColumns) {
            startUpdater();
        }

        isUpdating = false;
//        }
    }

    private void goToImage(int index) {
        tableModel.setSelected(index);

        int coords[] = tableModel.getRowColForIndex(index);

        // Gets current selected cell bounds.
        Rectangle rect = table.getCellRect(coords[0], coords[1], true);

        // Ensures item is visible.
        Point pos = jsPanel.getViewport().getViewPosition();
        rect.translate(-pos.x, -pos.y);
        jsPanel.getViewport().scrollRectToVisible(rect);

        repaint();
    }

    private void avgImage() {
        ImagesWindowFactory.captureFrame(ImageOperations.mean(tableModel.getAllItems()));
    }

    private void stdDevImage() {
        ImagesWindowFactory.captureFrame(ImageOperations.std_deviation(tableModel.getAllItems()));
    }

    private void setNormalized(boolean normalize) {
        jtbNormalize.setSelected(normalize);

        tableModel.setNormalized(normalize);

        startUpdater();
    }

    private void setUseGeometry(boolean useGeometry) {
        if (tableModel.containsGeometryInfo()) {
            ((MDTableModel) tableModel).setUseGeometry(useGeometry);

            updateTable();
        }
    }

    private void openAsStack() {
        ImagesWindowFactory.openTableAsImagePlus(tableModel);
    }

    private void openAs3D() {
        ImagesWindowFactory.openTableAs3D(tableModel);
    }

    private void saveAsMetadata(boolean all) {
        // Sets path and filename automatically.
        String filename = tableModel.getFilename() != null ? tableModel.getFilename() : "";
        fc.setSelectedFile(new File(forceExtension(filename, ".xmd")));

        if (fc.showSaveDialog(this) != JFileChooser.CANCEL_OPTION) {
            boolean response = true;
            if (fc.getSelectedFile().exists()) {
                response = JOptionPane.showConfirmDialog(null,
                        LABELS.MESSAGE_OVERWRITE_FILE,
                        LABELS.MESSAGE_OVERWRITE_FILE_TITLE,
                        JOptionPane.OK_CANCEL_OPTION,
                        JOptionPane.QUESTION_MESSAGE) == JOptionPane.OK_OPTION;
            }

            if (response) {
                String path = fc.getSelectedFile().getAbsolutePath();
                if (tableModel.saveAsMetadata(path, all)) {
                    JOptionPane.showMessageDialog(this, LABELS.MESSAGE_FILE_SAVED + path,
                            LABELS.MESSAGE_FILE_SAVED_TITLE, JOptionPane.INFORMATION_MESSAGE);
                }
            }
        }
    }

    private void saveAsStack(boolean all) {
        // Sets path and filename automatically.
        String filename = tableModel.getFilename() != null ? tableModel.getFilename() : "";
        fc.setSelectedFile(new File(forceExtension(filename, ".stk")));

        if (fc.showSaveDialog(this) != JFileChooser.CANCEL_OPTION) {
            boolean response = true;
            if (fc.getSelectedFile().exists()) {
                response = JOptionPane.showConfirmDialog(null,
                        LABELS.MESSAGE_OVERWRITE_FILE,
                        LABELS.MESSAGE_OVERWRITE_FILE_TITLE,
                        JOptionPane.OK_CANCEL_OPTION,
                        JOptionPane.QUESTION_MESSAGE) == JOptionPane.OK_OPTION;
            }

            if (response) {
                String path = fc.getSelectedFile().getAbsolutePath();
                if (tableModel.saveAsStack(path, all)) {
                    JOptionPane.showMessageDialog(this, LABELS.MESSAGE_FILE_SAVED + path,
                            LABELS.MESSAGE_FILE_SAVED_TITLE, JOptionPane.INFORMATION_MESSAGE);
                }
            }
        }
    }

    static String forceExtension(String filename, String ext) {
        int dot = filename.lastIndexOf(".");
        return filename.substring(0, dot) + ext;
    }

    public void pca() {
        String filename = tableModel.getFilename();

        try {
            MetaData md = new MetaData(filename);
            ImageDouble image = new ImageDouble();

            md.getPCAbasis(image);

            ImagesWindowFactory.captureFrame(ImageConverter.convertToImagej(image, "PCA: " + filename));
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }
    }

    public void fsc() {
        String filename = tableModel.getFilename();

        try {
            MetaData md = new MetaData(filename);

            ImagesWindowFactory.openFSCWindow(filename);
        } catch (Exception ex) {
            DEBUG.printException(ex);
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
        jtbNormalize = new javax.swing.JToggleButton();
        jtbUseGeometry = new javax.swing.JToggleButton();
        jpCenter = new javax.swing.JPanel();
        jpDisplay = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        jsZoom = new javax.swing.JSpinner();
        jcbShowLabels = new javax.swing.JCheckBox();
        jcbMDLabels = new javax.swing.JComboBox();
        jcbSortByLabel = new javax.swing.JCheckBox();
        jLabel2 = new javax.swing.JLabel();
        jsGoToImage = new javax.swing.JSpinner();
        jsPanel = new javax.swing.JScrollPane();
        table = new javax.swing.JTable();
        jpStructure = new javax.swing.JPanel();
        jcbAutoAdjustColumns = new javax.swing.JCheckBox();
        jLabel3 = new javax.swing.JLabel();
        jsRows = new javax.swing.JSpinner();
        jLabel4 = new javax.swing.JLabel();
        jsColumns = new javax.swing.JSpinner();

        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowOpened(java.awt.event.WindowEvent evt) {
                formWindowOpened(evt);
            }
        });
        addComponentListener(new java.awt.event.ComponentAdapter() {
            public void componentResized(java.awt.event.ComponentEvent evt) {
                formComponentResized(evt);
            }
        });

        toolBar.setRollover(true);

        jtbNormalize.setText(LABELS.BUTTON_NORMALIZE);
        jtbNormalize.setFocusable(false);
        jtbNormalize.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jtbNormalize.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jtbNormalize.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jtbNormalizeActionPerformed(evt);
            }
        });
        toolBar.add(jtbNormalize);

        jtbUseGeometry.setText(LABELS.LABEL_USE_GEOMETRY);
        jtbUseGeometry.setFocusable(false);
        jtbUseGeometry.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jtbUseGeometry.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jtbUseGeometry.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jtbUseGeometryActionPerformed(evt);
            }
        });
        toolBar.add(jtbUseGeometry);

        getContentPane().add(toolBar, java.awt.BorderLayout.NORTH);

        jpCenter.setLayout(new java.awt.BorderLayout());

        jLabel1.setText(LABELS.LABEL_ZOOM);
        jpDisplay.add(jLabel1);

        jsZoom.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(1), Integer.valueOf(1), null, Integer.valueOf(1)));
        jsZoom.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jsZoomStateChanged(evt);
            }
        });
        jpDisplay.add(jsZoom);

        jcbShowLabels.setText(LABELS.LABEL_SHOW_LABELS);
        jcbShowLabels.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jcbShowLabelsActionPerformed(evt);
            }
        });
        jpDisplay.add(jcbShowLabels);

        jcbMDLabels.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jcbMDLabelsItemStateChanged(evt);
            }
        });
        jpDisplay.add(jcbMDLabels);

        jcbSortByLabel.setText(LABELS.LABEL_SORT_BY_LABEL);
        jcbSortByLabel.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jcbSortByLabelActionPerformed(evt);
            }
        });
        jpDisplay.add(jcbSortByLabel);

        jLabel2.setText(LABELS.LABEL_GO2IMAGE);
        jpDisplay.add(jLabel2);

        jsGoToImage.setValue(1);
        jsGoToImage.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jsGoToImageStateChanged(evt);
            }
        });
        jpDisplay.add(jsGoToImage);

        jpCenter.add(jpDisplay, java.awt.BorderLayout.NORTH);

        table.setModel(tableModel);
        table.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
        table.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                tableMouseClicked(evt);
            }
        });
        jsPanel.setViewportView(table);

        jpCenter.add(jsPanel, java.awt.BorderLayout.CENTER);

        jcbAutoAdjustColumns.setText(LABELS.LABEL_AUTO_AJUST_COLUMNS);
        jcbAutoAdjustColumns.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jcbAutoAdjustColumnsStateChanged(evt);
            }
        });
        jpStructure.add(jcbAutoAdjustColumns);

        jLabel3.setText(LABELS.LABEL_ROWS);
        jpStructure.add(jLabel3);

        jsRows.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jsRowsStateChanged(evt);
            }
        });
        jpStructure.add(jsRows);

        jLabel4.setText(LABELS.LABEL_COLUMNS);
        jpStructure.add(jLabel4);

        jsColumns.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jsColumnsStateChanged(evt);
            }
        });
        jpStructure.add(jsColumns);

        jpCenter.add(jpStructure, java.awt.BorderLayout.SOUTH);

        getContentPane().add(jpCenter, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jsZoomStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jsZoomStateChanged
        setZoom((Integer) jsZoom.getValue());
    }//GEN-LAST:event_jsZoomStateChanged
    private void jcbShowLabelsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jcbShowLabelsActionPerformed
        setShowLabels(jcbShowLabels.isSelected());
    }//GEN-LAST:event_jcbShowLabelsActionPerformed
    private void formWindowOpened(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowOpened
        pack();
        ImagesWindowFactory.setConvenientSize(this);

        setInitialValues();
    }//GEN-LAST:event_formWindowOpened
    private void jsRowsStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jsRowsStateChanged
        if (!isUpdating) {
            setRows((Integer) jsRows.getValue());
        }
}//GEN-LAST:event_jsRowsStateChanged
    private void jsColumnsStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jsColumnsStateChanged
        if (!isUpdating) {
            setColumns((Integer) jsColumns.getValue());
        }
}//GEN-LAST:event_jsColumnsStateChanged
    private void jsGoToImageStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jsGoToImageStateChanged
        goToImage((Integer) jsGoToImage.getValue() - 1);
}//GEN-LAST:event_jsGoToImageStateChanged
    private void formComponentResized(java.awt.event.ComponentEvent evt) {//GEN-FIRST:event_formComponentResized
        startUpdater();
    }//GEN-LAST:event_formComponentResized
    private void tableMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_tableMouseClicked
        int view_row = table.rowAtPoint(evt.getPoint());
        int view_col = table.columnAtPoint(evt.getPoint());

        if (evt.getButton() == MouseEvent.BUTTON1) {  // Left click.
            if (evt.getClickCount() > 1) {
                Object item = table.getValueAt(view_row, view_col);

                if (item instanceof AbstractTableImageItem) {
                    ImagesWindowFactory.captureFrame(((AbstractTableImageItem) item).getImagePlus());
                }
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

    private void jtbNormalizeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jtbNormalizeActionPerformed
        setNormalized(jtbNormalize.isSelected());
    }//GEN-LAST:event_jtbNormalizeActionPerformed

    private void jcbMDLabelsItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jcbMDLabelsItemStateChanged
        if (evt.getStateChange() == ItemEvent.DESELECTED) {
            updateLabelToShow(jcbMDLabels.getSelectedIndex());
        }
    }//GEN-LAST:event_jcbMDLabelsItemStateChanged

    private void jcbSortByLabelActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jcbSortByLabelActionPerformed
        tableModel.setSorting(jcbSortByLabel.isSelected());
        updateTable();
    }//GEN-LAST:event_jcbSortByLabelActionPerformed

    private void jtbUseGeometryActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jtbUseGeometryActionPerformed
        setUseGeometry(jtbUseGeometry.isSelected());
}//GEN-LAST:event_jtbUseGeometryActionPerformed

    private void jcbAutoAdjustColumnsStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jcbAutoAdjustColumnsStateChanged
        setAutoAdjustColumns(jcbAutoAdjustColumns.isSelected());
}//GEN-LAST:event_jcbAutoAdjustColumnsStateChanged
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JCheckBox jcbAutoAdjustColumns;
    private javax.swing.JComboBox jcbMDLabels;
    private javax.swing.JCheckBox jcbShowLabels;
    private javax.swing.JCheckBox jcbSortByLabel;
    private javax.swing.JPanel jpCenter;
    private javax.swing.JPanel jpDisplay;
    protected javax.swing.JPanel jpStructure;
    protected javax.swing.JSpinner jsColumns;
    protected javax.swing.JSpinner jsGoToImage;
    private javax.swing.JScrollPane jsPanel;
    protected javax.swing.JSpinner jsRows;
    protected javax.swing.JSpinner jsZoom;
    private javax.swing.JToggleButton jtbNormalize;
    private javax.swing.JToggleButton jtbUseGeometry;
    private javax.swing.JTable table;
    private javax.swing.JToolBar toolBar;
    // End of variables declaration//GEN-END:variables

    class TableUpdater extends TimerTask {

        @Override
        public void run() {
            updateTable();
        }
    }

    class JMenuBarTable extends JMenuBar {

        protected JMenu jmSave = new JMenu(LABELS.LABEL_TABLE_SAVE);
        protected JMenuItem jmiSaveAsMetadata = new JMenuItem(LABELS.LABEL_TABLE_SAVE_AS_METADATA);
        protected JMenuItem jmiSaveAsStack = new JMenuItem(LABELS.LABEL_TABLE_SAVE_AS_STACK);
        protected JMenuItem jmiSaveSelectionAsMetadata = new JMenuItem(LABELS.LABEL_TABLE_SAVE_SELECTION_AS_METADATA);
        protected JMenuItem jmiSaveSelectionAsStack = new JMenuItem(LABELS.LABEL_TABLE_SAVE_SELECTION_AS_STACK);
        protected JMenu jmStatistics = new JMenu(LABELS.LABEL_MENU_STATISTICS);
        protected JMenuItem jmiAVG = new JMenuItem(LABELS.BUTTON_MEAN);
        protected JMenuItem jmiSTDEV = new JMenuItem(LABELS.BUTTON_STD_DEVIATION);
        protected JMenuItem jmiPCA = new JMenuItem(LABELS.BUTTON_PCA);
        protected JMenuItem jmiFSC = new JMenuItem(LABELS.BUTTON_FSC);
        protected JMenu jmOpenAs = new JMenu(LABELS.LABEL_MENU_OPEN_AS);
        protected JMenuItem jmiOpenAsStack = new JMenuItem(LABELS.BUTTON_TO_STACK);
        protected JMenuItem jmiOpenAsMetadata = new JMenuItem(LABELS.BUTTON_OPEN_AS_TABLE);
        protected JMenuItem jmiOpenAs3D = new JMenuItem(LABELS.OPERATION_OPEN_AS_3D);

        public JMenuBarTable() {
            super();

            add(jmSave);
            jmSave.add(jmiSaveAsMetadata);
            jmSave.add(jmiSaveAsStack);
            jmSave.addSeparator();
            jmSave.add(jmiSaveSelectionAsMetadata);
            jmSave.add(jmiSaveSelectionAsStack);

            jmiSaveAsMetadata.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    saveAsMetadata(true);
                }
            });

            jmiSaveSelectionAsMetadata.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    saveAsMetadata(false);
                }
            });

            jmiSaveAsStack.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    saveAsStack(true);
                }
            });

            jmiSaveSelectionAsStack.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    saveAsStack(false);
                }
            });

            add(jmStatistics);
            jmStatistics.add(jmiAVG);
            jmStatistics.add(jmiSTDEV);
            jmStatistics.add(jmiPCA);
            jmStatistics.add(jmiFSC);

            jmiAVG.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    avgImage();
                }
            });

            jmiSTDEV.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    stdDevImage();
                }
            });

            jmiPCA.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    pca();
                }
            });

            jmiFSC.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    fsc();
                }
            });

            add(jmOpenAs);
            jmOpenAs.add(jmiOpenAs3D);
            jmOpenAs.add(jmiOpenAsMetadata);
            jmOpenAs.add(jmiOpenAsStack);

            jmiOpenAs3D.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    openAs3D();
                }
            });

            jmiOpenAsMetadata.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    ImagesWindowFactory.openFileAsMetadata(tableModel.getFilename());
                }
            });

            jmiOpenAsStack.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    openAsStack();
                }
            });
        }
    }

    class JPopUpMenuTable extends JPopupMenu {

        protected Point location;
        protected JMenuItem jmiEnable = new JMenuItem(LABELS.LABEL_TABLE_ENABLE);
        protected JMenuItem jmiDisable = new JMenuItem(LABELS.LABEL_TABLE_DISABLE);
        protected JMenuItem jmiEnableAll = new JMenuItem(LABELS.LABEL_TABLE_ENABLE_ALL);
        protected JMenuItem jmiDisableAll = new JMenuItem(LABELS.LABEL_TABLE_DISABLE_ALL);
        protected JMenuItem jmiEnableFrom = new JMenuItem(LABELS.LABEL_TABLE_ENABLE_FROM);
        protected JMenuItem jmiEnableTo = new JMenuItem(LABELS.LABEL_TABLE_ENABLE_TO);
        protected JMenuItem jmiDisableFrom = new JMenuItem(LABELS.LABEL_TABLE_DISABLE_FROM);
        protected JMenuItem jmiDisableTo = new JMenuItem(LABELS.LABEL_TABLE_DISABLE_TO);

        public JPopUpMenuTable() {
            add(jmiEnable);
            add(jmiDisable);
            add(new JSeparator());
            add(jmiEnableAll);
            add(jmiDisableAll);
            add(new JSeparator());
            add(jmiEnableFrom);
            add(jmiEnableTo);
            add(new JSeparator());
            add(jmiDisableFrom);
            add(jmiDisableTo);

            jmiEnable.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, InputEvent.CTRL_DOWN_MASK));
            jmiDisable.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_D, InputEvent.CTRL_DOWN_MASK));
            jmiEnableAll.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_A, InputEvent.CTRL_DOWN_MASK));
            jmiDisableAll.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N, InputEvent.CTRL_DOWN_MASK));

//            jmiEnableFrom.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F, InputEvent.SHIFT_DOWN_MASK));
//            jmiDisableFrom.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F, InputEvent.CTRL_DOWN_MASK));
//            jmiEnableTo.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T, InputEvent.SHIFT_DOWN_MASK));
//            jmiDisableTo.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T, InputEvent.CTRL_DOWN_MASK));

            table.addKeyListener(new KeyListener() {

                public void keyTyped(KeyEvent e) {
                }

                public void keyPressed(KeyEvent e) {
                }

                public void keyReleased(KeyEvent e) {
                    if (e.isControlDown()) {
                        if (e.getKeyCode() == KeyEvent.VK_E) {
                            tableModel.enableSelectedItems();
                        }
                        if (e.getKeyCode() == KeyEvent.VK_D) {
                            tableModel.disableSelectedItems();
                        }
                        if (e.getKeyCode() == KeyEvent.VK_A) {
                            tableModel.enableAllItems();
                        }
                        if (e.getKeyCode() == KeyEvent.VK_N) {
                            tableModel.disableAllItems();
                        }
                    }
                }
            });

            jmiEnable.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    tableModel.enableSelectedItems();
                }
            });

            jmiDisable.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    tableModel.disableSelectedItems();
                }
            });

            jmiEnableAll.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    tableModel.enableAllItems();
                }
            });

            jmiDisableAll.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    tableModel.disableAllItems();
                }
            });

            jmiEnableFrom.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    int row = table.rowAtPoint(location);
                    int col = table.columnAtPoint(location);

                    tableModel.setEnabledFrom(row, col, true);
                }
            });

            jmiEnableTo.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    int row = table.rowAtPoint(location);
                    int col = table.columnAtPoint(location);

                    tableModel.setEnabledTo(row, col, true);
                }
            });

            jmiDisableFrom.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    int row = table.rowAtPoint(location);
                    int col = table.columnAtPoint(location);

                    tableModel.setEnabledFrom(row, col, false);
                }
            });

            jmiDisableTo.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    int row = table.rowAtPoint(location);
                    int col = table.columnAtPoint(location);

                    tableModel.setEnabledTo(row, col, false);
                }
            });
        }

        public void show(Component cmpnt, Point location) {
            this.location = location;
            show(cmpnt, location.x, location.y);
        }
    }
}
