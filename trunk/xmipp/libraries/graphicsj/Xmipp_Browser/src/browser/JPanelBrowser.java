/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */

/*
 * JPanelBrowser.java
 *
 * Created on 13-ene-2010, 16:02:05
 */
package browser;

import browser.windows.ImagesWindowFactory;
import browser.files.FileBrowser;
import browser.files.JComboBoxRootsModel;
import browser.files.FileListRenderer;
import browser.imageitems.listitems.FileItem;
import browser.files.FilterFilesModel;
import browser.imageitems.listitems.AbstractImageItem;
import browser.imageitems.listitems.SelFileItem;
import browser.imageitems.listitems.XmippImageItem;
import browser.windows.ImageWindowOperations;
import browser.windows.StackWindowOperations;
import ij.IJ;
import java.awt.BorderLayout;
import java.awt.event.KeyEvent;
import java.io.File;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import java.util.Vector;

/**
 *
 * @author Juanjo Vega
 */
public class JPanelBrowser extends JPanel {

    // Creates icons manager so icons will be ready to be used.
    protected boolean SHOW_PREVIEWS = true;
    protected FilterFilesModel filteringFilesModel;
    protected FileListRenderer fileListRenderer = new FileListRenderer();
    protected JComboBoxRootsModel comboBoxRootsModel = new JComboBoxRootsModel();
    protected JSearchBox searchBox = new JSearchBox(LABELS.LABEL_FILTER);

    /** Creates new form JPanelBrowser */
    public JPanelBrowser(String work_dir) {
        initComponents();

        filteringFilesModel = new FilterFilesModel(work_dir);
        filteringFilesModel.setFilteringLabel(jlFiltering);
        jlFilterFiles.setModel(filteringFilesModel);
        jlFilterFiles.installJTextField(searchBox.getTextField());
        jlFilterFiles.setCellRenderer(fileListRenderer);

        updateTitle();

        jcbRoots.setSelectedItem(filteringFilesModel.getCurrentRoot());

        jpFileBrowser.add(searchBox, BorderLayout.SOUTH);

        //jpImageInfo.setVisible(false);
    }

    protected void updateListStatus() {
        jlFilterFiles.setSelectedIndex(0);
        jlFilterFiles.ensureIndexIsVisible(0);
        updateTitle();
        updatePreview();
    }

    protected void updateTitle() {
        ((TitledBorder) jspFilesFiltering.getBorder()).setTitle(filteringFilesModel.getCurrentDirectory());
        jspFilesFiltering.updateUI();
    }

    // Updates preview with current item.
    public void updatePreview() {
        if (SHOW_PREVIEWS) {
            if (jlFilterFiles.getSelectedIndex() >= 0) {  // To avoid exceptions...
                Object item = jlFilterFiles.getSelectedValue();
                if (item instanceof AbstractImageItem) {
                    AbstractImageItem imageItem = (AbstractImageItem) item;

                    ImagePlus preview = imageItem.getPreview(ICONS_MANAGER.PREVIEW_WIDTH, ICONS_MANAGER.PREVIEW_HEIGHT);
                    if (preview != null) {
                        jpImageInfo.setPreview(preview.getImage(), imageItem.getImageInfo());
                        preview.close();
                    } else {
                        jpImageInfo.clearPreview();
                    }
                }
            }
        }

        // Shows / Hide preview panel.
        jpImageInfo.setVisible(SHOW_PREVIEWS);
    }

    // When files are opened usgin the "open" button.
    protected void openFileList(int indices[]) {
        for (int i = 0; i < indices.length; i++) {
            openFile(indices[i], indices.length == 1);   // Avoids going into directories for multiple selections.
        }
    }

    // When files are opened by clicking them or pressing the "enter" key.
    protected void openClickedFileList(int indices[]) {
        for (int i = 0; i < indices.length; i++) {
            openFile(indices[i], false);
        }
    }

    protected void openFile(int index, boolean openDir) {
        FileItem item = (FileItem) filteringFilesModel.getElementAt(index);
        File file = item.getFile();

        if (file.isDirectory()) {
            if (openDir) {
                filteringFilesModel.changeDirectory(index);
                updateListStatus();
            }
        } else {
            if (FileBrowser.canOpenFile(file)) {
                if (item instanceof SelFileItem) {
                    ImagesWindowFactory.openImage((SelFileItem) item);
                } else if (item instanceof XmippImageItem) {
                    XmippImageItem xmippItem = (XmippImageItem) item;
                    int n = 0;

                    if (xmippItem.dimension.nimages > 0) {
                        n = getImageIndexFromUser(xmippItem);
                    }

                    if (n >= 0) {
                        ImagesWindowFactory.openImage((XmippImageItem) item, n);
                    }
                } else {
                    IJ.error("@TODO Open with IJ (Xmipp can't open this file).");
                }
            } else {
                JOptionPane.showMessageDialog(this,
                        LABELS.MESSAGE_MEMORY_ERROR(file.length(), IJ.maxMemory()),
                        LABELS.TITLE_ERROR, JOptionPane.ERROR_MESSAGE);
            }
        }
    }

    // If there are more than une image, asks user for the one to open.
    private static int getImageIndexFromUser(XmippImageItem item) {
        int image = -1;

        if (item.dimension.nimages > 1) {
            final GenericDialog dialog = new GenericDialog("Select image from: " + item.getFileName());
            final String[] indexes = new String[(int) item.dimension.nimages];
            for (int i = 0; i < item.dimension.nimages; i++) {
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
    }

    private void send2Table(Object items[]) {
        ImagesWindowFactory.openTable(items);
    }

    private void captureFrames() {
        int ids[] = WindowManager.getIDList();
        Vector<String> windows = new Vector<String>();

        if (ids != null) {
            for (int i = 0; i < ids.length; i++) {
                ImageWindow window = WindowManager.getImage(ids[i]).getWindow();

                if (!(window instanceof ImageWindowOperations) && !(window instanceof StackWindowOperations)) {
                    windows.add(window.getTitle());
                }
            }
        }

        GenericDialog gd = new GenericDialog("Capture frame");
        String titles[] = new String[windows.size()];
        windows.toArray(titles);
        gd.addChoice("Windows:", titles, windows.size() > 0 ? titles[0] : "");
        gd.showDialog();
        if (windows.size() > 0 && !gd.wasCanceled()) {
            String selected = gd.getNextChoice();
            //ImagesWindowFactory.openImageWindow(WindowManager.getImage(selected), false);
            IJ.error("@TODO Capture frame");
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

        jpFileBrowser = new javax.swing.JPanel();
        jcbRoots = new javax.swing.JComboBox();
        jPanel1 = new javax.swing.JPanel();
        jspFilesFiltering = new javax.swing.JScrollPane();
        jlFilterFiles = new browser.files.JListFilterFiles();
        jlFiltering = new javax.swing.JLabel();
        jpPreview = new javax.swing.JPanel();
        jpImageInfo = new browser.JPanelImageInfo();
        jcbPreview = new javax.swing.JCheckBox();
        jpButtons = new javax.swing.JPanel();
        jbOpen = new javax.swing.JButton();
        jbSend2Table = new javax.swing.JButton();
        jToolBar = new javax.swing.JToolBar();
        jbParent = new javax.swing.JButton();
        jbRefresh = new javax.swing.JButton();
        jbCaptureWindow = new javax.swing.JButton();

        setLayout(new java.awt.BorderLayout());

        jpFileBrowser.setLayout(new java.awt.BorderLayout());

        jcbRoots.setModel(comboBoxRootsModel);
        jcbRoots.setSelectedIndex(0);
        jcbRoots.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jcbRootsActionPerformed(evt);
            }
        });
        jpFileBrowser.add(jcbRoots, java.awt.BorderLayout.NORTH);

        jPanel1.setLayout(new java.awt.BorderLayout());

        jspFilesFiltering.setBorder(javax.swing.BorderFactory.createTitledBorder(""));

        jlFilterFiles.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jlFilterFilesMouseClicked(evt);
            }
        });
        jlFilterFiles.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                jlFilterFilesKeyReleased(evt);
            }
        });
        jspFilesFiltering.setViewportView(jlFilterFiles);

        jPanel1.add(jspFilesFiltering, java.awt.BorderLayout.CENTER);
        jPanel1.add(jlFiltering, java.awt.BorderLayout.SOUTH);

        jpFileBrowser.add(jPanel1, java.awt.BorderLayout.CENTER);

        add(jpFileBrowser, java.awt.BorderLayout.CENTER);

        jpPreview.setLayout(new java.awt.BorderLayout());
        jpPreview.add(jpImageInfo, java.awt.BorderLayout.CENTER);

        jcbPreview.setSelected(true);
        jcbPreview.setText(LABELS.LABEL_LOAD_PREVIEW);
        jcbPreview.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jcbPreviewActionPerformed(evt);
            }
        });
        jpPreview.add(jcbPreview, java.awt.BorderLayout.SOUTH);

        add(jpPreview, java.awt.BorderLayout.EAST);

        jbOpen.setText(LABELS.BUTTON_OPEN);
        jbOpen.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbOpenActionPerformed(evt);
            }
        });
        jpButtons.add(jbOpen);

        jbSend2Table.setText(LABELS.BUTTON_SEND2TABLE);
        jbSend2Table.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbSend2TableActionPerformed(evt);
            }
        });
        jpButtons.add(jbSend2Table);

        add(jpButtons, java.awt.BorderLayout.SOUTH);

        jToolBar.setFloatable(false);
        jToolBar.setRollover(true);

        jbParent.setIcon(ICONS_MANAGER.PARENT_DIRECTORY);
        jbParent.setText(LABELS.BUTTON_PARENT_DIRECTORY);
        jbParent.setFocusable(false);
        jbParent.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbParent.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbParent.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbParentActionPerformed(evt);
            }
        });
        jToolBar.add(jbParent);

        jbRefresh.setIcon(ICONS_MANAGER.REFRESH_DIRECTORY);
        jbRefresh.setText(LABELS.BUTTON_REFRESH_DIRECTORY);
        jbRefresh.setFocusable(false);
        jbRefresh.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbRefresh.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbRefresh.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbRefreshActionPerformed(evt);
            }
        });
        jToolBar.add(jbRefresh);

        jbCaptureWindow.setIcon(ICONS_MANAGER.CAPTURE_WINDOW);
        jbCaptureWindow.setText(LABELS.BUTTON_CAPTURE_WINDOW);
        jbCaptureWindow.setFocusable(false);
        jbCaptureWindow.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbCaptureWindow.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbCaptureWindow.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbCaptureWindowActionPerformed(evt);
            }
        });
        jToolBar.add(jbCaptureWindow);

        add(jToolBar, java.awt.BorderLayout.NORTH);
    }// </editor-fold>//GEN-END:initComponents

    private void jbRefreshActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbRefreshActionPerformed
        filteringFilesModel.refreshCurrentDirectory();
        updateListStatus();
    }//GEN-LAST:event_jbRefreshActionPerformed

    private void jbParentActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbParentActionPerformed
        filteringFilesModel.goParent();
        updateListStatus();
    }//GEN-LAST:event_jbParentActionPerformed

    private void jcbRootsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jcbRootsActionPerformed
        // Changes unit
        File file = (File) jcbRoots.getSelectedItem();

        if (filteringFilesModel.getCurrentRoot().compareTo(file) != 0) {
            filteringFilesModel.changeDirectory((File) jcbRoots.getSelectedItem());
            updateListStatus();
        }
    }//GEN-LAST:event_jcbRootsActionPerformed

    private void jbOpenActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbOpenActionPerformed
        openFileList(jlFilterFiles.getSelectedIndices());
    }//GEN-LAST:event_jbOpenActionPerformed

    private void jlFilterFilesKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_jlFilterFilesKeyReleased
        if (evt.getKeyCode() == KeyEvent.VK_ENTER) {
            openClickedFileList(jlFilterFiles.getSelectedIndices());
        }

        updatePreview();
    }//GEN-LAST:event_jlFilterFilesKeyReleased

    private void jlFilterFilesMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jlFilterFilesMouseClicked
        if (evt.getClickCount() == 2) {
            int index = jlFilterFiles.locationToIndex(evt.getPoint());
            jlFilterFiles.ensureIndexIsVisible(index);
            openClickedFileList(jlFilterFiles.getSelectedIndices());
        }

        updatePreview();
    }//GEN-LAST:event_jlFilterFilesMouseClicked

    private void jcbPreviewActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jcbPreviewActionPerformed
        SHOW_PREVIEWS = jcbPreview.isSelected();
        updatePreview();
    }//GEN-LAST:event_jcbPreviewActionPerformed

    private void jbSend2TableActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbSend2TableActionPerformed
        send2Table(jlFilterFiles.getSelectedValues());
    }//GEN-LAST:event_jbSend2TableActionPerformed

    private void jbCaptureWindowActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbCaptureWindowActionPerformed
        captureFrames();
    }//GEN-LAST:event_jbCaptureWindowActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel jPanel1;
    private javax.swing.JToolBar jToolBar;
    private javax.swing.JButton jbCaptureWindow;
    private javax.swing.JButton jbOpen;
    private javax.swing.JButton jbParent;
    private javax.swing.JButton jbRefresh;
    private javax.swing.JButton jbSend2Table;
    private javax.swing.JCheckBox jcbPreview;
    private javax.swing.JComboBox jcbRoots;
    private browser.files.JListFilterFiles jlFilterFiles;
    private javax.swing.JLabel jlFiltering;
    private javax.swing.JPanel jpButtons;
    private javax.swing.JPanel jpFileBrowser;
    private browser.JPanelImageInfo jpImageInfo;
    private javax.swing.JPanel jpPreview;
    private javax.swing.JScrollPane jspFilesFiltering;
    // End of variables declaration//GEN-END:variables
}
