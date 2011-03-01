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
import browser.imageitems.listitems.SelFileItem;
import browser.imageitems.listitems.XmippImageItem;
import browser.table.normalization.iNormalizeListener;
import browser.utils.TaskTimer;
import ij.IJ;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import javax.swing.JFrame;
import javax.swing.JSpinner;
import javax.swing.JSpinner.DefaultEditor;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameVolumeTable extends JFrame implements iNormalizeListener {

    protected TaskTimer zoomTimer;    // Timer for zoom.
    private boolean isUpdating;
    protected int initialRows, initialColumns;

    /** Creates new form JFrameImagesTable */
    public JFrameVolumeTable() {
        this(0, 0);
    }

/*    private void startAbortTask() {
        TimerTask abort = new Abort();
        Timer timer = new Timer();
        timer.scheduleAtFixedRate(abort, 10000, 3000);
        System.err.println(" *** Started thread to avoid hangs!!!...");
    }*/

    public JFrameVolumeTable(int initialRows, int initialColumns) {
        super();

//        startAbortTask();

        this.initialRows = initialRows;
        this.initialColumns = initialColumns;

        initComponents();

        ((JSpinner.NumberEditor) jsRows.getEditor()).getTextField().setColumns(7);
        ((JSpinner.NumberEditor) jsColumns.getEditor()).getTextField().setColumns(7);
        ((JSpinner.NumberEditor) jsZoom.getEditor()).getTextField().setColumns(7);
        ((JSpinner.NumberEditor) jsGoToImage.getEditor()).getTextField().setColumns(7);

        jPanelTable.setShowLabels(false);   // No labels for volumes

        jsZoom.setFocusable(true);
        jsRows.setFocusable(true);
        jsColumns.setFocusable(true);
        jsGoToImage.setFocusable(true);

        ((DefaultEditor) jsRows.getEditor()).getTextField().addKeyListener(new KeyListener() {

            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                if (e.getKeyCode() == KeyEvent.VK_ENTER
                        || e.getKeyCode() == KeyEvent.VK_UP
                        || e.getKeyCode() == KeyEvent.VK_DOWN) {
                    setRows();
                }
            }
        });

        ((DefaultEditor) jsColumns.getEditor()).getTextField().addKeyListener(new KeyListener() {

            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                if (e.getKeyCode() == KeyEvent.VK_ENTER
                        || e.getKeyCode() == KeyEvent.VK_UP
                        || e.getKeyCode() == KeyEvent.VK_DOWN) {
                    setColumns();
                }
            }
        });

        ((DefaultEditor) jsGoToImage.getEditor()).getTextField().addKeyListener(new KeyListener() {

            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                if (e.getKeyCode() == KeyEvent.VK_ENTER
                        || e.getKeyCode() == KeyEvent.VK_UP
                        || e.getKeyCode() == KeyEvent.VK_DOWN) {
                    goToImage();
                }
            }
        });

        zoomTimer = new TaskTimer(new Runnable() {

            public void run() {
                jPanelTable.setZoom((Integer) jsZoom.getValue());
            }
        });
    }

    protected void setInitialDisplayValues() {
        // Calculates and sets the best initial zoom.
        if (jPanelTable.getItemsCount() > 0) {
            int zoom = 100;

            Dimension cellSize = new Dimension(jPanelTable.getCellWidth(), jPanelTable.getCellHeight());
            Dimension CELL_SIZE = jPanelTable.getCellSize();
            Rectangle PANEL_SIZE = jPanelTable.getVisibleRect();

            int bigger_panel_dimension, cell_dimension;

            if (PANEL_SIZE.width > PANEL_SIZE.height) {
                bigger_panel_dimension = PANEL_SIZE.width;
                cell_dimension = CELL_SIZE.width;
            } else {
                bigger_panel_dimension = PANEL_SIZE.height;
                cell_dimension = CELL_SIZE.height;
            }

            double scale = (double) bigger_panel_dimension / (double) cell_dimension;

            if (scale < 1) {
                zoom *= scale;
            }

            setZoom(zoom);

            // Calculates and sets columns.
            int items = PANEL_SIZE.width / (cellSize.width > 0 ? cellSize.width : 1);
            setColumns(items > jPanelTable.getItemsCount() ? jPanelTable.getItemsCount() : items);
        } else {
            jsRows.setEnabled(false);
            jsColumns.setEnabled(false);
            jsGoToImage.setEnabled(false);
        }
    }

    public void addImageItem(XmippImageItem itemImage) {
        jPanelTable.addImageItem(itemImage);

        setTitle(LABELS.TITLE_TABLE_WINDOW(itemImage));
    }

    /*    public void addImageItem(XmippImageItem itemImage) {
    for (int i = 0; i < itemImage.dimension.nimages; i++) {
    jPanelTable.addImageItem(itemImage, i);
    }

    setTitle(LABELS.TITLE_TABLE_WINDOW(itemImage));
    }*/
//
//    public void addImageItem(XmippImageItem itemImage, int n) {
//        jPanelTable.addImageItem(itemImage, n);
//
//        setTitle(LABELS.TITLE_TABLE_WINDOW(itemImage));
//    }
    public void addImageItem(SelFileItem itemImage) {
        jPanelTable.addImageItem(itemImage);

        setTitle(LABELS.TITLE_TABLE_WINDOW(itemImage));
    }
//    public void addVolume(XmippImageItem itemImage) {
//        jPanelTable.addImageItem(itemImage, 0);
//        setTitle(LABELS.TITLE_TABLE_WINDOW_VOLUME(itemImage));
//    }
//
//    public void addSelFile(SelFileItem itemSel) {
//        jPanelTable.addSelFile(itemSel);
//        setTitle(
//                LABELS.TITLE_TABLE_WINDOW_SELFILE(
//                itemSel));
//    }

    public int getImagesCount() {
        return jPanelTable.getItemsCount();
    }

    protected int getZoom() {
        return (Integer) jsZoom.getValue();
    }

    protected void setRows(int rows) {
        jsRows.setValue(rows);
        setRows();
    }

    protected void setColumns(int columns) {
        jsColumns.setValue(columns);
        setColumns();
    }

    protected void setRows() {
        int rows = (Integer) jsRows.getValue();

        jPanelTable.setRows(rows);

        jsColumns.setValue(jPanelTable.getColumnsCount());
    }

    protected void setColumns() {
        int cols = (Integer) jsColumns.getValue();

        jPanelTable.setColumns(cols);

        jsRows.setValue(jPanelTable.getRowsCount());
    }

    protected void goToImage() {
        jPanelTable.setImageSelected((Integer) jsGoToImage.getValue() - 1);
    }

    public void normalize(double min, double max) {
        jPanelTable.normalize(min, max);
    }

    public void setNormalizedAuto() {
        jPanelTable.setNormalizedAuto();
    }

    public void disableNormalization() {
        jPanelTable.disableNormalization();
    }

    private void setZoom(int zoom) {
        isUpdating = true;
        jsZoom.setValue(zoom);
        isUpdating = false;

        // Starts zoom zoomTimer.
        zoomTimer.start();
    }

    private void send2Stack() {
//        ImagesWindowFactory.openTableAsStack(jPanelTable.getItems());
        IJ.error("@TODO SEnd to stack");
    }

    private void autoAdjustColumns() {
        jPanelTable.autoAdjustColumns();
    }

    public void setAutoAdjustColumns(boolean autoAdjustColumns) {
        jPanelTable.setAutoAdjustColumns(autoAdjustColumns);
    }

    public void setTableSize(int rows, int columns) {
        jPanelTable.setTableSize(rows, columns);

        isUpdating = true;  // Avoids upodating events ;)
        jsRows.setValue(jPanelTable.getRowsCount());
        jsColumns.setValue(jPanelTable.getColumnsCount());
        isUpdating = false; // Restore listener processing.
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
        jbToStack = new javax.swing.JButton();
        jbSelectAll = new javax.swing.JButton();
        jbInvertSelection = new javax.swing.JButton();
        jbRefresh = new javax.swing.JButton();
        jpContent = new javax.swing.JPanel();
        jpZoom = new javax.swing.JPanel();
        jsZoom = new javax.swing.JSpinner();
        jcbAutoAdjustColumns = new javax.swing.JCheckBox();
        jPanelTable = new browser.table.JPanelTable();
        jpControls = new javax.swing.JPanel();
        jsRows = new javax.swing.JSpinner();
        jsColumns = new javax.swing.JSpinner();
        jsGoToImage = new javax.swing.JSpinner();

        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosed(java.awt.event.WindowEvent evt) {
                formWindowClosed(evt);
            }
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

        jbSelectAll.setText(LABELS.BUTTON_SELECT_ALL);
        jbSelectAll.setFocusable(false);
        jbSelectAll.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbSelectAll.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbSelectAll.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbSelectAllActionPerformed(evt);
            }
        });
        toolBar.add(jbSelectAll);

        jbInvertSelection.setText(LABELS.BUTTON_INVERT_SELECTION);
        jbInvertSelection.setFocusable(false);
        jbInvertSelection.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbInvertSelection.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbInvertSelection.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbInvertSelectionActionPerformed(evt);
            }
        });
        toolBar.add(jbInvertSelection);

        jbRefresh.setText(LABELS.BUTTON_REFRESH);
        jbRefresh.setFocusable(false);
        jbRefresh.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbRefresh.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbRefresh.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbRefreshActionPerformed(evt);
            }
        });
        toolBar.add(jbRefresh);

        getContentPane().add(toolBar, java.awt.BorderLayout.NORTH);

        jpContent.setLayout(new java.awt.BorderLayout());

        jpZoom.setLayout(new java.awt.BorderLayout());

        jsZoom.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(0), Integer.valueOf(0), null, Integer.valueOf(1)));
        jsZoom.setBorder(javax.swing.BorderFactory.createTitledBorder(LABELS.LABEL_ZOOM));
        jsZoom.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jsZoomStateChanged(evt);
            }
        });
        jpZoom.add(jsZoom, java.awt.BorderLayout.WEST);

        jcbAutoAdjustColumns.setText(LABELS.LABEL_AUTO_AJUST_COLUMNS);
        jcbAutoAdjustColumns.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jcbAutoAdjustColumnsActionPerformed(evt);
            }
        });
        jpZoom.add(jcbAutoAdjustColumns, java.awt.BorderLayout.CENTER);

        jpContent.add(jpZoom, java.awt.BorderLayout.NORTH);
        jpContent.add(jPanelTable, java.awt.BorderLayout.CENTER);

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

        jpContent.add(jpControls, java.awt.BorderLayout.SOUTH);

        getContentPane().add(jpContent, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void formWindowClosed(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosed
        jPanelTable.clear();
    }//GEN-LAST:event_formWindowClosed
    private void jsGoToImageStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jsGoToImageStateChanged
        Integer value = (Integer) jsGoToImage.getValue();

        if (value < 1) {
            value = 1;
        } else if (value > jPanelTable.getItemsCount()) {
            value = jPanelTable.getItemsCount();
        }

        jsGoToImage.setValue(value);

        goToImage();
    }//GEN-LAST:event_jsGoToImageStateChanged
    private void jbSelectAllActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbSelectAllActionPerformed
        jPanelTable.selectAll();
    }//GEN-LAST:event_jbSelectAllActionPerformed
    private void jbInvertSelectionActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbInvertSelectionActionPerformed
        jPanelTable.invertSelection();
    }//GEN-LAST:event_jbInvertSelectionActionPerformed
    private void jbRefreshActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbRefreshActionPerformed
        jPanelTable.refreshData();
    }//GEN-LAST:event_jbRefreshActionPerformed
    private void jsZoomStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jsZoomStateChanged
        setZoom((Integer) jsZoom.getValue());
    }//GEN-LAST:event_jsZoomStateChanged
    private void formWindowOpened(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowOpened
        setInitialDisplayValues();

        if (initialRows > 0) {
            if (initialColumns > 0) {
                setTableSize(initialRows, initialColumns);
            } else {
                setRows(initialRows);
            }
        } else if (initialColumns > 0) {
            setColumns(initialColumns);
        }
    }//GEN-LAST:event_formWindowOpened
    private void formComponentResized(java.awt.event.ComponentEvent evt) {//GEN-FIRST:event_formComponentResized
        jPanelTable.refreshCacheSize();
        autoAdjustColumns();
    }//GEN-LAST:event_formComponentResized
    private void jsRowsStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jsRowsStateChanged
        if (!isUpdating) {
            isUpdating = true;

            int rows = (Integer) jsRows.getValue();

            if (rows < 1) {
                rows = 1;
            } else if (rows > jPanelTable.getItemsCount()) {
                rows = jPanelTable.getItemsCount();
            }

            setRows(rows);

            isUpdating = false;
        }
    }//GEN-LAST:event_jsRowsStateChanged
    private void jsColumnsStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jsColumnsStateChanged
        if (!isUpdating) {
            isUpdating = true;

            int columns = (Integer) jsColumns.getValue();

            if (columns < 1) {
                columns = 1;
            } else if (columns > jPanelTable.getItemsCount()) {
                columns = jPanelTable.getItemsCount();
            }

            setColumns(columns);

            isUpdating = false;
        }
    }//GEN-LAST:event_jsColumnsStateChanged

    private void jbToStackActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbToStackActionPerformed
        send2Stack();
    }//GEN-LAST:event_jbToStackActionPerformed

    private void jcbAutoAdjustColumnsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jcbAutoAdjustColumnsActionPerformed
        setAutoAdjustColumns(jcbAutoAdjustColumns.isSelected());
    }//GEN-LAST:event_jcbAutoAdjustColumnsActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    protected browser.table.JPanelTable jPanelTable;
    protected javax.swing.JButton jbInvertSelection;
    protected javax.swing.JButton jbRefresh;
    protected javax.swing.JButton jbSelectAll;
    protected javax.swing.JButton jbToStack;
    private javax.swing.JCheckBox jcbAutoAdjustColumns;
    private javax.swing.JPanel jpContent;
    protected javax.swing.JPanel jpControls;
    private javax.swing.JPanel jpZoom;
    protected javax.swing.JSpinner jsColumns;
    protected javax.swing.JSpinner jsGoToImage;
    protected javax.swing.JSpinner jsRows;
    protected javax.swing.JSpinner jsZoom;
    protected javax.swing.JToolBar toolBar;
    // End of variables declaration//GEN-END:variables
}
