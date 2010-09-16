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
import browser.imageitems.FileImageItem;
import browser.imageitems.SelFileItem;
import browser.table.normalization.iNormalizeListener;
import browser.utils.TaskTimer;
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
    protected boolean isUpdating;

    /** Creates new form JFrameImagesTable */
    public JFrameVolumeTable() {
        super();

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
        if (jPanelTable.getImagesCount() > 0) {
            int zoom = 100;

            Dimension cellSize = /*jPanelTable.getCellSize();*/ new Dimension(jPanelTable.getCellWidth(), jPanelTable.getCellHeight());
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
            setColumns(items > jPanelTable.getImagesCount() ? jPanelTable.getImagesCount() : items);
        } else {
            jsRows.setEnabled(false);
            jsColumns.setEnabled(false);
            jsGoToImage.setEnabled(false);
        }
    }

    public void addVolume(FileImageItem itemImage) {
        jPanelTable.addVolume(itemImage);

        setTitle(LABELS.TITLE_TABLE_WINDOW_VOLUME(itemImage));
    }

    public void addVolume(SelFileItem itemImage) {
        jPanelTable.addVolume(itemImage);

        setTitle(LABELS.TITLE_TABLE_WINDOW_SELFILE(itemImage));
    }

    public int getImagesCount() {
        return jPanelTable.getImagesCount();
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

        jsColumns.setValue(jPanelTable.getColumns());
    }

    protected void setColumns() {
        int cols = (Integer) jsColumns.getValue();

        jPanelTable.setColumns(cols);

        jsRows.setValue(jPanelTable.getRows());
    }

    protected void goToImage() {
        jPanelTable.setImageSelected((Integer) jsGoToImage.getValue());
        jPanelTable.repaint();
    }

    public void normalize(double min, double max) {
        jPanelTable.normalize(min, max);
    }

    public void normalizeAuto() {
        jPanelTable.normalizeAuto();
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

    /*    public static void main(String args[]) {
    int H = 403, W = 453;
    int h = 816, w = 816;

    int zoom = 100;

    double scale = (double) H / (double) h;
    //double scale = (double) h / (double) H;
    if (scale < 1) {
    zoom *= scale;
    }

    System.out.println(" *** SCALE: " + scale);
    System.out.println(" *** ZOOM: " + zoom);
    }*/
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        toolBar = new javax.swing.JToolBar();
        jbSelectAll = new javax.swing.JButton();
        jbInvertSelection = new javax.swing.JButton();
        jbRefresh = new javax.swing.JButton();
        jpContent = new javax.swing.JPanel();
        jpZoom = new javax.swing.JPanel();
        jsZoom = new javax.swing.JSpinner();
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

        jsZoom.setModel(new javax.swing.SpinnerNumberModel(0, 0, 200, 1));
        jsZoom.setBorder(javax.swing.BorderFactory.createTitledBorder(LABELS.LABEL_ZOOM));
        jsZoom.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jsZoomStateChanged(evt);
            }
        });
        jpZoom.add(jsZoom, java.awt.BorderLayout.WEST);

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



        if (value < 0) {
            value = 0;


        } else if (value > jPanelTable.getImagesCount() - 1) {
            value = jPanelTable.getImagesCount() - 1;


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
        /*        jSliderZoom.setValue((Integer) jsZoom.getValue());

        if (!isUpdating) {  // If it is updating the zoom, waits for zoomTimer to do it.
        setZoom();
        }*/
        setZoom((Integer) jsZoom.getValue());
    }//GEN-LAST:event_jsZoomStateChanged
    private void formWindowOpened(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowOpened
        setInitialDisplayValues();
    }//GEN-LAST:event_formWindowOpened
    private void formComponentResized(java.awt.event.ComponentEvent evt) {//GEN-FIRST:event_formComponentResized
        jPanelTable.refreshCacheSize();
    }//GEN-LAST:event_formComponentResized
    private void jsRowsStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jsRowsStateChanged
        int rows = (Integer) jsRows.getValue();

        if (rows < 1) {
            rows = 1;
        } else if (rows > jPanelTable.getImagesCount()) {
            rows = jPanelTable.getImagesCount();
        }

        setRows(rows);
    }//GEN-LAST:event_jsRowsStateChanged
    private void jsColumnsStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jsColumnsStateChanged
        int columns = (Integer) jsColumns.getValue();

        if (columns < 1) {
            columns = 1;
        } else if (columns > jPanelTable.getImagesCount()) {
            columns = jPanelTable.getImagesCount();
        }

        setColumns(columns);
    }//GEN-LAST:event_jsColumnsStateChanged
    // Variables declaration - do not modify//GEN-BEGIN:variables
    protected browser.table.JPanelTable jPanelTable;
    protected javax.swing.JButton jbInvertSelection;
    protected javax.swing.JButton jbRefresh;
    protected javax.swing.JButton jbSelectAll;
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
