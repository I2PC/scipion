/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
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
import browser.imageitems.FileImageItem;
import browser.imageitems.SelFileItem;
import ij.IJ;
import java.awt.BorderLayout;
import java.awt.event.KeyEvent;
import java.io.File;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import browser.table.JFrameImagesTable;
import browser.table.JFrameVolumeTable;
import ij.ImagePlus;
import xmipp.Sel_Reader;

/**
 *
 * @author Juanjo Vega
 */
public class JPanelBrowser extends JPanel {

    // Creates icons manager so icons will be ready to be used.
    protected ICONS_MANAGER ICONS_MANAGER = new ICONS_MANAGER();
    protected FilterFilesModel filteringFilesModel = new FilterFilesModel();
    protected FileListRenderer fileListRenderer = new FileListRenderer();
    protected JComboBoxRootsModel comboBoxRootsModel = new JComboBoxRootsModel();
    protected JSearchBox searchBox = new JSearchBox(LABELS.LABEL_FILTER);

    /** Creates new form JPanelBrowser */
    public JPanelBrowser() {
        initComponents();

        filteringFilesModel.setFilteringLabel(jlFiltering);
        jlFilterFiles.setModel(filteringFilesModel);
        jlFilterFiles.installJTextField(searchBox.getTextField());
        jlFilterFiles.setCellRenderer(fileListRenderer);

        updateTitle();

        jcbRoots.setSelectedItem(filteringFilesModel.getCurrentRoot());

        jpFileBrowser.add(searchBox, BorderLayout.SOUTH);

        jpImageInfo.setVisible(false);
    }

    protected void updateListStatus() {
        jlFilterFiles.setSelectedIndex(0);
        jlFilterFiles.ensureIndexIsVisible(0);
        updateTitle();
    }

    protected void updateTitle() {
        ((TitledBorder) jspFilesFiltering.getBorder()).setTitle(filteringFilesModel.getCurrentDirectory());
        jspFilesFiltering.updateUI();
    }

    // Updates preview with current item.
    public void updatePreview() {
        if (Model.SHOW_PREVIEWS) {
            if (jlFilterFiles.getSelectedIndex() >= 0) {  // To avoid exceptions...
                FileItem item = (FileItem) jlFilterFiles.getSelectedValue();

                ImagePlus preview = item.getPreview(ICONS_MANAGER.PREVIEW_WIDTH, ICONS_MANAGER.PREVIEW_HEIGHT);
                if (preview != null) {
                    jpImageInfo.setPreview(preview.getImage(), item.getFileInfo());
                    preview.close();
                }
            }
        }

        // Shows / Hide preview panel.
        jpImageInfo.setVisible(Model.SHOW_PREVIEWS);
    }

    protected void openFileList(int indices[]) {
        for (int i = 0; i < indices.length; i++) {
            openFile(indices[i], indices.length == 1);   // Avoids going into directories for multiple selections.
        }
    }

    protected void openFile(int index, boolean openDir) {
        File file = ((FileItem) filteringFilesModel.getElementAt(index)).getFile();

        if (file.isDirectory()) {
            if (openDir) {
                filteringFilesModel.changeDirectory(index);
                updateListStatus();
            }
        } else {
            if (FileBrowser.canOpenFile(file)) {
                ImagesWindowFactory.openImageWindow(file);
            } else {
                JOptionPane.showMessageDialog(this,
                        LABELS.MESSAGE_MEMORY_ERROR(file.length(), IJ.maxMemory()),
                        LABELS.TITLE_ERROR, JOptionPane.ERROR_MESSAGE);
            }
        }
    }

    private void send2Table(int indices[]) {
        JFrameImagesTable imagesTable = new JFrameImagesTable();

        for (int i = 0; i < indices.length; i++) {
            FileItem item = (FileItem) jlFilterFiles.getModel().getElementAt(indices[i]);
            if (item instanceof FileImageItem) {
                FileImageItem imageItem = (FileImageItem) item;
                if (imageItem.nslices > 1) {    // Volumes are opened in a independent table for each one
                    JFrameVolumeTable volumeTable = new JFrameVolumeTable();

                    volumeTable.addVolume(imageItem);

                    volumeTable.setVisible(true);
                } else {    // Images will be at the same table.
                    imagesTable.addImage(imageItem);
                }
            } else if (item instanceof SelFileItem) {
                JFrameVolumeTable volumeTable = new JFrameVolumeTable();

                volumeTable.addVolume((SelFileItem) item);

                volumeTable.setVisible(true);
            }
        }

        // If there was any image file.
        if (imagesTable.getImagesCount() > 0) {
            imagesTable.setVisible(true);
        }
    }

    private void send2Table_(int indices[]) {
        JFrameImagesTable imagesTable = new JFrameImagesTable();
        JFrameVolumeTable volumeTable = new JFrameVolumeTable();

        for (int i = 0; i < indices.length; i++) {
            FileItem item = (FileItem) jlFilterFiles.getModel().getElementAt(indices[i]);
            if (item instanceof FileImageItem) {
                FileImageItem imageItem = (FileImageItem) item;
                if (imageItem.nslices > 1) {    // Volumes are opened in a independent table for each one
                    volumeTable.addVolume(imageItem);

                    volumeTable.setVisible(true);
                } else {    // Images will be at the same table.
                    imagesTable.addImage(imageItem);
                }
            }
        }

        // If there was any image file.
        if (imagesTable.getImagesCount() > 0) {
            imagesTable.setVisible(true);
        }

        if (volumeTable.getImagesCount() > 0) {
            volumeTable.setVisible(true);
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
            openFileList(jlFilterFiles.getSelectedIndices());
        }

        updatePreview();
    }//GEN-LAST:event_jlFilterFilesKeyReleased

    private void jlFilterFilesMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jlFilterFilesMouseClicked
        if (evt.getClickCount() == 2) {
            int index = jlFilterFiles.locationToIndex(evt.getPoint());
            jlFilterFiles.ensureIndexIsVisible(index);
            openFileList(jlFilterFiles.getSelectedIndices());
        }

        updatePreview();
    }//GEN-LAST:event_jlFilterFilesMouseClicked

    private void jcbPreviewActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jcbPreviewActionPerformed
        Model.SHOW_PREVIEWS = jcbPreview.isSelected();
        updatePreview();
    }//GEN-LAST:event_jcbPreviewActionPerformed

    private void jbSend2TableActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbSend2TableActionPerformed
        send2Table(jlFilterFiles.getSelectedIndices());
    }//GEN-LAST:event_jbSend2TableActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel jPanel1;
    private javax.swing.JToolBar jToolBar;
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
