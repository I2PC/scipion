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

import browser.DEBUG;
import browser.LABELS;
import browser.imageitems.TableImageItem;
import browser.table.renderers.ImageRenderer;
import browser.table.renderers.RowHeaderRenderer;
import browser.windows.ImagesWindowFactory;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.util.Timer;
import java.util.TimerTask;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.JSeparator;
import javax.swing.KeyStroke;
import javax.swing.LookAndFeel;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.Border;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameImagesTable extends JFrame {//implements TableModelListener {

    private final static int DELAY_TO_UPDATE = 500;
    private ImagesTableModel tableModel;
    private ImagesTableColumnModel columnModel;
    private ImagesRowHeaderModel rowHeaderModel;
    private JList rowHeader;
    private ImageRenderer renderer = new ImageRenderer();
    private Timer updateTimer = new Timer(true);    // Timer for zoom.
    private TableUpdater tableUpdaterTask;    // Associated task for zoom timer.
    private boolean isUpdating;
    private boolean autoAdjustColumns = false;
    private JPopUpMenuTable jpopUpMenuTable;
    private JMenuBarTable jMenuBarTable;

    /** Creates new form JFrameImagesTable */
    public JFrameImagesTable(String filename) {
        this(filename, 0, 0);
    }

    public JFrameImagesTable(String filenames[]) {
        this(filenames, 0, 0);
    }

    // @TODO Clean initialRows and cols: It might not be used anymore.
    public JFrameImagesTable(String filename, int initialRows, int initialColumns) {
        super();

        tableModel = new ImagesTableModel(filename);
        postInit();
    }

    public JFrameImagesTable(String filenames[], int initialRows, int initialColumns) {
        this(filenames, null, initialRows, initialColumns);
    }

    public JFrameImagesTable(String filenames[], boolean enabled[]) {
        this(filenames, enabled, 0, 0);
    }

    public JFrameImagesTable(String filenames[], boolean enabled[], int initialRows, int initialColumns) {
        super();

        tableModel = new ImagesTableModel(filenames, enabled);
        postInit();
    }

    private void postInit() {
        jMenuBarTable = new JMenuBarTable();
        setJMenuBar(jMenuBarTable);

        columnModel = new ImagesTableColumnModel();

        initComponents();

        table.setColumnModel(columnModel);
        table.setDefaultRenderer(TableImageItem.class, renderer);

        setTitle(tableModel.getTitle());

        jpopUpMenuTable = new JPopUpMenuTable();

        // Sets limits for spinners.
        jsRows.setModel(new SpinnerNumberModel(1, 1, tableModel.getSize(), 1));
        jsColumns.setModel(new SpinnerNumberModel(1, 1, tableModel.getSize(), 1));
        jsGoToImage.setModel(new SpinnerNumberModel(0, 0, tableModel.getSize() - 1, 1));

        setRowHeader();

        // Stacks will be "auto-normalized".
        if (tableModel.isStack()) {
            setNormalized(true);
        }
    }

    private void setRowHeader() {
        rowHeaderModel = new ImagesRowHeaderModel(table, 1);

        rowHeader = new JList();
        rowHeader.setModel(rowHeaderModel);

        LookAndFeel.installColorsAndFont(rowHeader, "TableHeader.background",
                "TableHeader.foreground", "TableHeader.font");

        rowHeader.setCellRenderer(new RowHeaderRenderer());

        jsPanel.setRowHeaderView(rowHeader);
    }

    private void updateTable() {
//        if (table.isShowing()) {
        isUpdating = true;
        DEBUG.printMessage(" *** Updating table: " + System.currentTimeMillis());

        if (tableModel.getSize() > 0) {
            Dimension dimension = getCellSize();
            renderer.setPreferredSize(dimension);

            // Adjusts rows size.
            table.setRowHeight(dimension.height > 0 ? dimension.height : 1);
            rowHeader.setFixedCellHeight(table.getRowHeight());
//            rowHeader.setFixedCellWidth(rowHeader.getModel().getSize() * 2);
            columnModel.setWidth(dimension.width > 0 ? dimension.width : 1);

            autoAdjustColumns();    // If auto adjust columns is enabled, refresh!

            rowHeader.repaint();
            table.revalidate();
        }

        isUpdating = false;
//        }
    }

    private Dimension getCellSize() {
        int font_height = 0;
        if (renderer.isShowingLabels()) {
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
        renderer.setShowLabels(show);

//        startUpdater();
        updateTable();
    }

    private void autoAdjustColumns() {
        if (autoAdjustColumns) {
            tableModel.autoAdjustColumns(jsPanel.getViewportBorderBounds().width - jsPanel.getVerticalScrollBar().getWidth(),
                    table.getIntercellSpacing().width);

            jsRows.setValue(tableModel.getRowCount());
            jsColumns.setValue(tableModel.getColumnCount());
        }
    }

    private void setInitialValues() {
        if (tableModel.getSize() > 0) {
            setZoom((int) (getInitialScale() * 100));

            setAutoAdjustColumns(true);
        } else {
            jsZoom.setEnabled(false);
            jcbAutoAdjustColumns.setEnabled(false);
            jcbShowLabels.setEnabled(false);
            jsRows.setEnabled(false);
            jsColumns.setEnabled(false);
            jsGoToImage.setEnabled(false);
        }
    }

    private double getInitialScale() {
        double W = jsPanel.getVisibleRect().width;
        double w = ((TableImageItem) tableModel.getValueAt(0, 0)).getWidth();

        double scale = W / w;

        return scale > 1.0 ? 1.0 : scale;
    }

    private void setZoom(int zoom) {
        isUpdating = true;
        tableModel.setZoomScale(zoom / 100.0);

        jsZoom.setValue(zoom);

        startUpdater();
        isUpdating = false;
    }

    private void setRows(int rows) {
        isUpdating = true;
        tableModel.setRows(rows);

        jsColumns.setValue(tableModel.getColumnCount());

        startUpdater();
        isUpdating = false;
    }

    private void setColumns(int columns) {
        isUpdating = true;
        tableModel.setColumns(columns);

        jsRows.setValue(tableModel.getRowCount());

        startUpdater();
        isUpdating = false;
    }

    private void setAutoAdjustColumns(boolean autoAdjustColumns) {
        isUpdating = true;
        jcbAutoAdjustColumns.setSelected(autoAdjustColumns);

        this.autoAdjustColumns = autoAdjustColumns;

        if (autoAdjustColumns) {
            autoAdjustColumns();
        }

        jsColumns.setEnabled(!jcbAutoAdjustColumns.isSelected());
        jsRows.setEnabled(!jcbAutoAdjustColumns.isSelected());

        startUpdater();
        //updateTable();
        isUpdating = false;
    }

    private void goToImage(int index) {
        System.out.println(" @TODO: going to image: " + index);

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

    private void send2stack() {
        ImagesWindowFactory.openTableAsImagePlus(tableModel);
    }

    private void openAs3D() {
        ImagesWindowFactory.openTableAs3D(tableModel);
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {
        java.awt.GridBagConstraints gridBagConstraints;

        toolBar = new javax.swing.JToolBar();
        jtbNormalize = new javax.swing.JToggleButton();
        jbMean = new javax.swing.JButton();
        jbStdDev = new javax.swing.JButton();
        jbToStack = new javax.swing.JButton();
        jbTo3D = new javax.swing.JButton();
        jpCenter = new javax.swing.JPanel();
        jpZoom = new javax.swing.JPanel();
        jsZoom = new javax.swing.JSpinner();
        jcbAutoAdjustColumns = new javax.swing.JCheckBox();
        jcbShowLabels = new javax.swing.JCheckBox();
        jsPanel = new javax.swing.JScrollPane();
        table = new javax.swing.JTable();
        jpControls = new javax.swing.JPanel();
        jsRows = new javax.swing.JSpinner();
        jsColumns = new javax.swing.JSpinner();
        jsGoToImage = new javax.swing.JSpinner();

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

        jbMean.setText(LABELS.BUTTON_MEAN);
        jbMean.setFocusable(false);
        jbMean.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbMean.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbMean.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbMeanActionPerformed(evt);
            }
        });
        toolBar.add(jbMean);

        jbStdDev.setText(LABELS.BUTTON_STD_DEVIATION);
        jbStdDev.setFocusable(false);
        jbStdDev.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbStdDev.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbStdDev.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbStdDevActionPerformed(evt);
            }
        });
        toolBar.add(jbStdDev);

        jbToStack.setText(LABELS.BUTTON_TO_STACK);
        jbToStack.setFocusable(false);
        jbToStack.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbToStack.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbToStack.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbToStackActionPerformed(evt);
            }
        });
        toolBar.add(jbToStack);

        jbTo3D.setText(LABELS.OPERATION_OPEN_AS_3D);
        jbTo3D.setFocusable(false);
        jbTo3D.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbTo3D.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbTo3D.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbTo3DActionPerformed(evt);
            }
        });
        toolBar.add(jbTo3D);

        getContentPane().add(toolBar, java.awt.BorderLayout.NORTH);

        jpCenter.setLayout(new java.awt.BorderLayout());

        jpZoom.setLayout(new java.awt.GridBagLayout());

        jsZoom.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(1), Integer.valueOf(1), null, Integer.valueOf(1)));
        jsZoom.setBorder(javax.swing.BorderFactory.createTitledBorder(LABELS.LABEL_ZOOM));
        jsZoom.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jsZoomStateChanged(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.ipadx = 100;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
        jpZoom.add(jsZoom, gridBagConstraints);

        jcbAutoAdjustColumns.setText(LABELS.LABEL_AUTO_AJUST_COLUMNS);
        jcbAutoAdjustColumns.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jcbAutoAdjustColumnsActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.anchor = java.awt.GridBagConstraints.WEST;
        jpZoom.add(jcbAutoAdjustColumns, gridBagConstraints);

        jcbShowLabels.setText(LABELS.LABEL_SHOW_LABELS);
        jcbShowLabels.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jcbShowLabelsActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.anchor = java.awt.GridBagConstraints.WEST;
        jpZoom.add(jcbShowLabels, gridBagConstraints);

        jpCenter.add(jpZoom, java.awt.BorderLayout.NORTH);

        table.setModel(tableModel);
        table.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
        table.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                tableMouseClicked(evt);
            }
        });
        jsPanel.setViewportView(table);

        jpCenter.add(jsPanel, java.awt.BorderLayout.CENTER);

        jpControls.setLayout(new java.awt.GridLayout(1, 0));

        jsRows.setBorder(javax.swing.BorderFactory.createTitledBorder(LABELS.LABEL_ROWS));
        jsRows.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jsRowsStateChanged(evt);
            }
        });
        jpControls.add(jsRows);

        jsColumns.setBorder(javax.swing.BorderFactory.createTitledBorder(LABELS.LABEL_COLUMNS));
        jsColumns.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jsColumnsStateChanged(evt);
            }
        });
        jpControls.add(jsColumns);

        jsGoToImage.setBorder(javax.swing.BorderFactory.createTitledBorder(LABELS.LABEL_GO2IMAGE));
        jsGoToImage.setValue(1);
        jsGoToImage.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jsGoToImageStateChanged(evt);
            }
        });
        jpControls.add(jsGoToImage);

        jpCenter.add(jpControls, java.awt.BorderLayout.SOUTH);

        getContentPane().add(jpCenter, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jsZoomStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jsZoomStateChanged
        setZoom((Integer) jsZoom.getValue());
    }//GEN-LAST:event_jsZoomStateChanged
    private void jcbShowLabelsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jcbShowLabelsActionPerformed
        setShowLabels(jcbShowLabels.isSelected());
    }//GEN-LAST:event_jcbShowLabelsActionPerformed
    private void jcbAutoAdjustColumnsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jcbAutoAdjustColumnsActionPerformed
        setAutoAdjustColumns(jcbAutoAdjustColumns.isSelected());
    }//GEN-LAST:event_jcbAutoAdjustColumnsActionPerformed
    private void formWindowOpened(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowOpened
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
        goToImage((Integer) jsGoToImage.getValue());
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

                if (item instanceof TableImageItem) {
                    ImagesWindowFactory.captureFrame(((TableImageItem) item).getImagePlus());
                }
            } else {
                // Ctrl adds items to selection, otherwise previous ones are removed.
                if (!evt.isControlDown()) {
                    tableModel.clearSelection();
                }

                // Right click selects an item (but doesn't deselect it)
                tableModel.toggleSelected(view_row, view_col);

                table.repaint();
            }
        } else if (evt.getButton() == MouseEvent.BUTTON3) {  // Right click.
            table.setRowSelectionInterval(view_row, view_row);
            table.setColumnSelectionInterval(view_col, view_col);

            if (!evt.isControlDown()) {
                tableModel.clearSelection();
            }

            tableModel.setSelected(view_row, view_col, true);

            jpopUpMenuTable.show(evt.getComponent(), evt.getX(), evt.getY());
        }
    }//GEN-LAST:event_tableMouseClicked

    private void jbMeanActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbMeanActionPerformed
        avgImage();
    }//GEN-LAST:event_jbMeanActionPerformed

    private void jbStdDevActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbStdDevActionPerformed
        stdDevImage();
    }//GEN-LAST:event_jbStdDevActionPerformed

    private void jtbNormalizeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jtbNormalizeActionPerformed
        setNormalized(jtbNormalize.isSelected());
    }//GEN-LAST:event_jtbNormalizeActionPerformed

    private void jbTo3DActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbTo3DActionPerformed
        openAs3D();
    }//GEN-LAST:event_jbTo3DActionPerformed

    private void jbToStackActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbToStackActionPerformed
        send2stack();
    }//GEN-LAST:event_jbToStackActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton jbMean;
    private javax.swing.JButton jbStdDev;
    private javax.swing.JButton jbTo3D;
    private javax.swing.JButton jbToStack;
    private javax.swing.JCheckBox jcbAutoAdjustColumns;
    private javax.swing.JCheckBox jcbShowLabels;
    private javax.swing.JPanel jpCenter;
    protected javax.swing.JPanel jpControls;
    private javax.swing.JPanel jpZoom;
    protected javax.swing.JSpinner jsColumns;
    protected javax.swing.JSpinner jsGoToImage;
    private javax.swing.JScrollPane jsPanel;
    protected javax.swing.JSpinner jsRows;
    protected javax.swing.JSpinner jsZoom;
    private javax.swing.JToggleButton jtbNormalize;
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

        protected JMenu jmiSave = new JMenu(LABELS.LABEL_TABLE_SAVE);
        protected JMenuItem jmiSaveAsImages = new JMenuItem(LABELS.LABEL_TABLE_SAVE_AS_IMAGES);
        protected JMenuItem jmiSaveAsStack = new JMenuItem(LABELS.LABEL_TABLE_SAVE_AS_STACK);
        protected JMenuItem jmiSaveAsSelfile = new JMenuItem(LABELS.LABEL_TABLE_SAVE_AS_SELFILE);

        public JMenuBarTable() {
            super();

            add(jmiSave);
            jmiSave.add(jmiSaveAsImages);
            jmiSave.add(jmiSaveAsStack);
            jmiSave.add(jmiSaveAsSelfile);

            jmiSaveAsImages.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    /*                    if (!saveAsImages()) {
                    IJ.error(LABELS.MESSAGE_NO_ITEMS_SELECTED);
                    //                    JOptionPane.showMessageDialog(getParent(), LABELS.MESSAGE_NO_ITEMS_SELECTED);
                    }*/
                    System.out.println("@TODO Save as IMAGES file");
                }
            });

            jmiSaveAsStack.addActionListener(
                    new ActionListener() {

                        public void actionPerformed(ActionEvent e) {
                            /*                    if (!saveAsStack()) {
                            JOptionPane.showMessageDialog(getParent(), LABELS.MESSAGE_NO_ITEMS_SELECTED);
                            }*/
                            System.out.println("@TODO Save as STACK file");
                        }
                    });

            jmiSaveAsSelfile.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    /*                    if (!saveAsSelFile()) {
                    JOptionPane.showMessageDialog(getParent(), LABELS.MESSAGE_NO_ITEMS_SELECTED);
                    }*/
                    System.out.println("@TODO Save as SEL file");
                }
            });
        }
    }

    class JPopUpMenuTable extends JPopupMenu {

        protected JMenuItem jmiEnable = new JMenuItem(LABELS.LABEL_TABLE_ENABLE);
        protected JMenuItem jmiDisable = new JMenuItem(LABELS.LABEL_TABLE_DISABLE);
        protected JMenuItem jmiEnableAll = new JMenuItem(LABELS.LABEL_TABLE_ENABLE_ALL);
        protected JMenuItem jmiDisableAll = new JMenuItem(LABELS.LABEL_TABLE_DISABLE_ALL);

        public JPopUpMenuTable() {
            add(jmiEnable);
            add(jmiDisable);
            add(new JSeparator());
            add(jmiEnableAll);
            add(jmiDisableAll);
            add(new JSeparator());

            jmiEnable.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, InputEvent.CTRL_DOWN_MASK));
            jmiDisable.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_D, InputEvent.CTRL_DOWN_MASK));
            jmiEnableAll.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_A, InputEvent.CTRL_DOWN_MASK));
            jmiDisableAll.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N, InputEvent.CTRL_DOWN_MASK));

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

        }
    }
}
