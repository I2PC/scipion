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

import browser.LABELS;
import browser.imageitems.TableImageItem;
import browser.table.renderers.ImageRenderer;
import browser.table.renderers.RowHeaderRenderer;
import java.awt.Dimension;
import java.awt.Insets;
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
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.JSeparator;
import javax.swing.KeyStroke;
import javax.swing.border.Border;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameImagesTable extends JFrame {//implements TableModelListener {

    private final static int ZOOMDELAY = 2000;
    private ImagesTableModel tableModel;
    private ImagesTableColumnModel columnModel;
//    private ImagesRowHeaderModel rowHeaderModel;
    private JList rowHeader;
    private RowHeaderRenderer rowHeaderRenderer;
    private ImageRenderer renderer = new ImageRenderer();
    private Timer updateTimer = new Timer(true);    // Timer for zoom.
    private TableUpdater tableUpdaterTask;    // Associated task for zoom timer.
    private boolean isUpdating;
    private boolean autoAdjustColumns = false;
    private JPopUpMenuTable jpopUpMenuTable;

    /** Creates new form JFrameImagesTable */
    public JFrameImagesTable(String filename) {
        this(filename, 0, 0);
    }

    public JFrameImagesTable(String filenames[]) {
        this(filenames, 0, 0);
    }

    public JFrameImagesTable(String filename, int initialRows, int initialColumns) {
        super();

        tableModel = new ImagesTableModel(filename);
        postInit();
    }

    public JFrameImagesTable(String filenames[], int initialRows, int initialColumns) {
        super();

        tableModel = new ImagesTableModel(filenames);
        postInit();
    }

    private void postInit() {
        columnModel = new ImagesTableColumnModel();

        initComponents();

        jtImages.setColumnModel(columnModel);
        jtImages.setDefaultRenderer(TableImageItem.class, renderer);

        setTitle(tableModel.getTitle());

        jpopUpMenuTable = new JPopUpMenuTable();

//        tableModel.addTableModelListener(this);
//        setRowHeader();
    }

/*    public void setRowHeader() {
        rowHeaderModel = new ImagesRowHeaderModel(tableModel);

        rowHeader = new JList();
        rowHeader.setModel(rowHeaderModel);

        LookAndFeel.installColorsAndFont(rowHeader, "TableHeader.background",
                "TableHeader.foreground", "TableHeader.font");

        rowHeader.setFixedCellHeight(jtImages.getRowHeight());

        rowHeaderRenderer = new RowHeaderRenderer();
        rowHeader.setCellRenderer(rowHeaderRenderer);

        //rowHeader.addListSelectionListener(this);

        jspTable.setRowHeaderView(rowHeader);
    }*/

    /*private void updateZoom() {
    tableModel.setZoomScale((Integer) jsZoom.getValue());
    }*/
    /*private void updateRows() {
    tableModel.setRows((Integer) jsRows.getValue());
    }

    private void updateColumns() {
    tableModel.setColumns((Integer) jsColumns.getValue());
    }*/
    /*
    public void tableChanged(TableModelEvent e) {
    System.out.println("Table Changed! > " + e.getType() + " TYPE_UP: " + TableModelEvent.UPDATE);
    //        if (!isUpdating) {
    startUpdater();
    //        }
    }*/
    private void updateTable() {
        isUpdating = true;
        System.out.println(" *** Updating table: " + System.currentTimeMillis());

        if (tableModel.getSize() > 0) {
            Dimension dimension = getCellSize();
            renderer.setPreferredSize(dimension);

            // Adjusts rows size.
            jtImages.setRowHeight(dimension.height > 0 ? dimension.height : 1);
            columnModel.setWidth(dimension.width > 0 ? dimension.width : 1);

            autoAdjustColumns();    // If auto adjust columns is enabled, refresh!

//            tableModel.fireTableStructureChanged();
            jtImages.revalidate();
        }

        isUpdating = false;
    }

    private Dimension getCellSize() {
        int font_height = 0;
        if (renderer.isShowingLabels()) {
            font_height = renderer.getFontMetrics(renderer.getFont()).getHeight();
            font_height += renderer.getIconTextGap();  // Adds the extra gap.
            font_height -= jtImages.getIntercellSpacing().height;   // Removes spacing.
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
        updateTimer.schedule(tableUpdaterTask, ZOOMDELAY);
    }

    private void setShowLabels(boolean show) {
        renderer.setShowLabels(show);

        startUpdater();
    }

    private void autoAdjustColumns() {
        if (autoAdjustColumns) {
            tableModel.autoAdjustColumns(jspTable.getVisibleRect().width);

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
        double W = jspTable.getVisibleRect().width;
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
        isUpdating = false;
    }

    private void goToImage() {
        System.out.println(" @TODO: going to image...");
    }

    private void avgImage() {
        ImageOperations.average(tableModel.getAllItems()).show();
    }

    private void stdDevImage() {
        ImageOperations.std_deviation(tableModel.getAllItems()).show();
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
        jbAverage = new javax.swing.JButton();
        jbStdDev = new javax.swing.JButton();
        jpCenter = new javax.swing.JPanel();
        jpZoom = new javax.swing.JPanel();
        jsZoom = new javax.swing.JSpinner();
        jcbAutoAdjustColumns = new javax.swing.JCheckBox();
        jcbShowLabels = new javax.swing.JCheckBox();
        jspTable = new javax.swing.JScrollPane();
        jtImages = new javax.swing.JTable();
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

        jbAverage.setText(LABELS.BUTTON_AVERAGE);
        jbAverage.setFocusable(false);
        jbAverage.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbAverage.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbAverage.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbAverageActionPerformed(evt);
            }
        });
        toolBar.add(jbAverage);

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

        getContentPane().add(toolBar, java.awt.BorderLayout.NORTH);

        jpCenter.setLayout(new java.awt.BorderLayout());

        jpZoom.setLayout(new java.awt.GridBagLayout());

        jsZoom.setBorder(javax.swing.BorderFactory.createTitledBorder(LABELS.LABEL_ZOOM));
        jsZoom.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jsZoomStateChanged(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.ipadx = 100;
        jpZoom.add(jsZoom, gridBagConstraints);

        jcbAutoAdjustColumns.setText(LABELS.LABEL_AUTO_AJUST_COLUMNS);
        jcbAutoAdjustColumns.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jcbAutoAdjustColumnsActionPerformed(evt);
            }
        });
        jpZoom.add(jcbAutoAdjustColumns, new java.awt.GridBagConstraints());

        jcbShowLabels.setText(LABELS.LABEL_SHOW_LABELS);
        jcbShowLabels.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jcbShowLabelsActionPerformed(evt);
            }
        });
        jpZoom.add(jcbShowLabels, new java.awt.GridBagConstraints());

        jpCenter.add(jpZoom, java.awt.BorderLayout.NORTH);

        jtImages.setModel(tableModel);
        jtImages.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
        jtImages.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jtImagesMouseClicked(evt);
            }
        });
        jspTable.setViewportView(jtImages);

        jpCenter.add(jspTable, java.awt.BorderLayout.CENTER);

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
        Integer value = (Integer) jsGoToImage.getValue();

        if (value < 1) {
            value = 1;
        } else if (value > tableModel.getSize()) {
            value = tableModel.getSize();
        }

        jsGoToImage.setValue(value);

        goToImage();
}//GEN-LAST:event_jsGoToImageStateChanged
    private void formComponentResized(java.awt.event.ComponentEvent evt) {//GEN-FIRST:event_formComponentResized
        startUpdater();
        //updateTable();
    }//GEN-LAST:event_formComponentResized
    private void jtImagesMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jtImagesMouseClicked
        if (evt.getButton() == MouseEvent.BUTTON1) {  // Left click.
            if (evt.getClickCount() > 1) {
                //openSelection();
                System.out.println(" >>> Open selection.");
            } else {
                int col = jtImages.columnAtPoint(evt.getPoint());
                int row = jtImages.rowAtPoint(evt.getPoint());

                // Ctrl adds items to selection, otherwise previous ones are removed.
                if (!evt.isControlDown()) {
                    tableModel.clearSelection();
                }

                // Right click selects an item (but doesn't deselect it)
                tableModel.toggleSelected(row, col);

                jtImages.repaint();
            }
        } else if (evt.getButton() == MouseEvent.BUTTON3) {  // Right click.
            if (tableModel.getSelectedItems() != null) {
                jpopUpMenuTable.show(evt.getComponent(), evt.getX(), evt.getY());
            }
        }
    }//GEN-LAST:event_jtImagesMouseClicked

    private void jbAverageActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbAverageActionPerformed
        avgImage();
    }//GEN-LAST:event_jbAverageActionPerformed

    private void jbStdDevActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbStdDevActionPerformed
        stdDevImage();
    }//GEN-LAST:event_jbStdDevActionPerformed

    private void jtbNormalizeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jtbNormalizeActionPerformed
        if (jtbNormalize.isSelected()) {
            System.out.println("@ENABLING normalization");
            tableModel.setNormalized();
        } else {
            System.out.println("@DISABLING normalization");
            tableModel.disableNormalization();
        }

        tableModel.printStuff();

        startUpdater();
    }//GEN-LAST:event_jtbNormalizeActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton jbAverage;
    private javax.swing.JButton jbStdDev;
    private javax.swing.JCheckBox jcbAutoAdjustColumns;
    private javax.swing.JCheckBox jcbShowLabels;
    private javax.swing.JPanel jpCenter;
    protected javax.swing.JPanel jpControls;
    private javax.swing.JPanel jpZoom;
    protected javax.swing.JSpinner jsColumns;
    protected javax.swing.JSpinner jsGoToImage;
    protected javax.swing.JSpinner jsRows;
    protected javax.swing.JSpinner jsZoom;
    private javax.swing.JScrollPane jspTable;
    private javax.swing.JTable jtImages;
    private javax.swing.JToggleButton jtbNormalize;
    private javax.swing.JToolBar toolBar;
    // End of variables declaration//GEN-END:variables

    class TableUpdater extends TimerTask {

        @Override
        public void run() {
            updateTable();
        }
    }

    class JPopUpMenuTable extends JPopupMenu {

        protected JMenuItem jmiEnable = new JMenuItem(LABELS.LABEL_TABLE_ENABLE);
        protected JMenuItem jmiDisable = new JMenuItem(LABELS.LABEL_TABLE_DISABLE);
        protected JMenuItem jmiEnableAll = new JMenuItem(LABELS.LABEL_TABLE_ENABLE_ALL);
        protected JMenuItem jmiDisableAll = new JMenuItem(LABELS.LABEL_TABLE_DISABLE_ALL);
        protected JMenuItem jmiSaveAsImages = new JMenuItem(LABELS.LABEL_TABLE_SAVE_AS_IMAGES);
        protected JMenuItem jmiSaveAsStack = new JMenuItem(LABELS.LABEL_TABLE_SAVE_AS_STACK);
        protected JMenuItem jmiSaveAsSelfile = new JMenuItem(LABELS.LABEL_TABLE_SAVE_AS_SELFILE);

        public JPopUpMenuTable() {
            add(jmiEnable);
            add(jmiDisable);
            add(new JSeparator());
            add(jmiEnableAll);
            add(jmiDisableAll);
            add(new JSeparator());
            add(jmiSaveAsImages);
            add(jmiSaveAsStack);
            add(jmiSaveAsSelfile);

            jmiEnable.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, InputEvent.CTRL_DOWN_MASK));
            jmiDisable.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_D, InputEvent.CTRL_DOWN_MASK));
            jmiEnableAll.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_A, InputEvent.CTRL_DOWN_MASK));
            jmiDisableAll.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N, InputEvent.CTRL_DOWN_MASK));

            jtImages.addKeyListener(new KeyListener() {

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

            jmiSaveAsImages.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    /*                    if (!saveAsImages()) {
                    IJ.error(LABELS.MESSAGE_NO_ITEMS_SELECTED);
                    //                    JOptionPane.showMessageDialog(getParent(), LABELS.MESSAGE_NO_ITEMS_SELECTED);
                    }*/
                    System.out.println("@TODO Save as SEL file");
                }
            });

            jmiSaveAsStack.addActionListener(
                    new ActionListener() {

                        public void actionPerformed(ActionEvent e) {
                            /*                    if (!saveAsStack()) {
                            JOptionPane.showMessageDialog(getParent(), LABELS.MESSAGE_NO_ITEMS_SELECTED);
                            }*/
                            System.out.println("@TODO Save as SEL file");
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
}
