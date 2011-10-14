/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */

/*
 * JPanelBrowser.java
 *
 * Created on 13-ene-2010, 16:02:05
 */
package browser.filebrowsers;

import browser.InfiniteProgressPanel;
import browser.ICONS_MANAGER;
import browser.LABELS;
import browser.filebrowsers.model.FileBrowser;
import java.awt.BorderLayout;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import browser.filebrowsers.renderers.FileListRenderer;
import browser.filebrowsers.model.ListModelFilesBrowser;
import browser.imageitems.listitems.AbstractImageItem;
import browser.imageitems.listitems.FileItem;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import ij.ImagePlus;
import java.awt.event.KeyEvent;
import java.io.File;
import java.util.LinkedList;
import javax.swing.ListModel;
import javax.swing.ListSelectionModel;

/**
 *
 * @author Juanjo Vega
 */
public class JPanelXmippBrowser extends JPanel {

    // Creates icons manager so icons will be ready to be used.
    protected boolean SHOW_PREVIEWS = true;
    protected int previewWidth = ICONS_MANAGER.DEFAULT_PREVIEW_WIDTH;
    protected int previewHeight = ICONS_MANAGER.DEFAULT_PREVIEW_HEIGHT;
    protected ListModelFilesBrowser listModelFilesList;
    protected FileListRenderer fileListRenderer = new FileListRenderer();
    protected JSearchBox searchBox = new JSearchBox(LABELS.LABEL_FILTER);

    /** Creates new form JPanelBrowser */
    public JPanelXmippBrowser(String folder) {
        this(folder, "");
    }

    public JPanelXmippBrowser(final String folder, final String expression) {
        super();

        initComponents();

        Thread t = new Thread(new Runnable() {

            public void run() {
                listModelFilesList = new ListModelFilesBrowser(folder);
                listModelFilesList.setFilteringLabel(jlFiltering);
                jlFileFilter.setModel((ListModel) listModelFilesList);

                jlFileFilter.installJTextField(searchBox.getTextField());
                jlFileFilter.setCellRenderer(fileListRenderer);
                jpFileBrowser.add(searchBox, BorderLayout.SOUTH);

                // Sets initial expression.
                searchBox.textField.setText(expression);

                updatePath();
            }
        });

        t.start();
    }

    public void setSingleSelection(boolean single) {
        if (single) {
            jlFileFilter.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        } else {
            jlFileFilter.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        }
    }

    protected void updateListStatus() {
        jlFileFilter.setSelectedIndex(0);
        jlFileFilter.ensureIndexIsVisible(0);
        updatePath();
        updatePreview();
    }

    private void updatePath() {
        ((TitledBorder) jspFilesFiltering.getBorder()).setTitle(
                listModelFilesList.getCurrentDirectory());
        jspFilesFiltering.updateUI();
    }

    public void updatePreview() {
        clearPreview();

        if (SHOW_PREVIEWS) {
            final InfiniteProgressPanel progressPanel = new InfiniteProgressPanel("Calculating preview...");
            //final Component previousGlassPane = getRootPane().getGlassPane();
            getRootPane().setGlassPane(progressPanel);
            progressPanel.start();

            Thread t = new Thread(new Runnable() {

                public void run() {

                    if (jlFileFilter.getSelectedIndex() >= 0) {  // To avoid exceptions...
                        Object item = jlFileFilter.getSelectedValue();
                        if (item instanceof AbstractImageItem) {
                            setPreview((AbstractImageItem) item);
                        }
                    }

                    // Shows / Hide preview panel.
                    jpImageInfo.setVisible(SHOW_PREVIEWS);

                    progressPanel.setVisible(false);
                    progressPanel.stop();
                }
            });

            t.start();
        }
    }

    protected void setPreview(AbstractImageItem imageItem) {
        ImagePlus preview = null;

        if (imageItem != null) {
            preview = getPreview(imageItem);
        }

        if (preview != null) {
            jpImageInfo.setPreview(preview.getImage(), imageItem.getImageInfo());
        } else {
            jpImageInfo.clearPreview();
        }
    }

    public void setPreviewSize(int width, int height) {
        previewWidth = width;
        previewHeight = height;

        jpImageInfo.setPreviewSize(width, height);
    }

    protected ImagePlus getPreview(AbstractImageItem item) {
        return item.getPreview(previewWidth, previewHeight);
    }

    protected void clearPreview() {
        jpImageInfo.clearPreview();
    }

    public void refreshCurrentDirectory() {
        listModelFilesList.refreshCurrentDirectory();
        updateListStatus();
    }

    public void goParent() {
        listModelFilesList.goParent();
        updateListStatus();
    }

    public Object[] getSelectedValues() {
        return jlFileFilter.getSelectedValues();
    }

    // When a file is open by clicking it or pressing the "enter" key.
    protected void openSelectedFile() {
        if (jlFileFilter.getSelectedIndex() == 0) {
            goParent();
        } else {
            FileItem item = (FileItem) jlFileFilter.getSelectedValue();

            if (item.isDirectory()) {
                listModelFilesList.changeDirectory(item.getFile());
                updateListStatus();
            } else {
                openFileAsDefault(item);
            }
        }
    }

    public void openSelectedFilesAsImages() {
        openFilesAsImage(jlFileFilter.getSelectedValues());
    }

    public void openSelectedFilesAsTable() {
        openFilesAsTable(jlFileFilter.getSelectedValues());
    }

    protected void openFileAsDefault(FileItem item) {
        if (FileBrowser.hasEnoughMemory(item.getFile())) {
            ImagesWindowFactory.openFileAsDefault(item.getAbsoluteFileName());
        } else {
            IJ.showMessage(LABELS.TITLE_ERROR,
                    LABELS.MESSAGE_MEMORY_ERROR(item.getFile().length(), IJ.maxMemory()));
        }
    }

    private void openFilesAsImage(Object items[]) {
        for (int i = 0; i < items.length; i++) {
            openFileAsImage((FileItem) items[i]);
        }
    }

    private void openFileAsImage(FileItem item) {
        File file = item.getFile();

        if (!item.isDirectory()) {
            if (FileBrowser.hasEnoughMemory(file)) {
                ImagesWindowFactory.openFileAsImage(item.getAbsoluteFileName());
            } else {
                IJ.showMessage(LABELS.TITLE_ERROR,
                        LABELS.MESSAGE_MEMORY_ERROR(file.length(), IJ.maxMemory()));
            }
        }
    }
    /*
     * Separates single images from the rest of files.
     * Single images are opened in the same table together.
     * Stacks, volumes and metadata files are opened in indepedent tables.
     */

    private void openFilesAsTable(Object items[]) {
        LinkedList<String> images = new LinkedList<String>();

        for (int i = 0; i < items.length; i++) {
            Object object = items[i];
            if (object instanceof AbstractImageItem) {
                AbstractImageItem abstractItem = (AbstractImageItem) object;

                if (abstractItem.isSingleImage()) {
                    images.add(abstractItem.getAbsoluteFileName());
                } else {
                    ImagesWindowFactory.openFileAsGallery(abstractItem.getAbsoluteFileName());
                }
            } else {
                clearPreview();
            }
        }

        if (!images.isEmpty()) {
            String imagesArray[] = images.toArray(new String[images.size()]);
            ImagesWindowFactory.openFilesAsGallery(imagesArray, true);
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
        jpCenter = new javax.swing.JPanel();
        jspFilesFiltering = new javax.swing.JScrollPane();
        jlFileFilter = new browser.filebrowsers.JListFileFilter();
        jlFiltering = new javax.swing.JLabel();
        jpPreview = new javax.swing.JPanel();
        jpImageInfo = new browser.filebrowsers.JPanelImageInfo();
        jcbPreview = new javax.swing.JCheckBox();

        setLayout(new java.awt.BorderLayout());

        jpFileBrowser.setLayout(new java.awt.BorderLayout());

        jpCenter.setLayout(new java.awt.BorderLayout());

        jspFilesFiltering.setBorder(javax.swing.BorderFactory.createTitledBorder(""));

        jlFileFilter.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jlFileFilterMouseClicked(evt);
            }
        });
        jlFileFilter.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                jlFileFilterKeyReleased(evt);
            }
        });
        jspFilesFiltering.setViewportView(jlFileFilter);

        jpCenter.add(jspFilesFiltering, java.awt.BorderLayout.CENTER);
        jpCenter.add(jlFiltering, java.awt.BorderLayout.SOUTH);

        jpFileBrowser.add(jpCenter, java.awt.BorderLayout.CENTER);

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
    }// </editor-fold>//GEN-END:initComponents

    private void jlFileFilterKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_jlFileFilterKeyReleased
        if (evt.getKeyCode() == KeyEvent.VK_ENTER) {
            openSelectedFile();
        }

        updatePreview();
    }//GEN-LAST:event_jlFileFilterKeyReleased

    private void jlFileFilterMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jlFileFilterMouseClicked
        if (evt.getClickCount() == 2) {
            openSelectedFile();
        }

        updatePreview();
    }//GEN-LAST:event_jlFileFilterMouseClicked

    private void jcbPreviewActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jcbPreviewActionPerformed
        SHOW_PREVIEWS = jcbPreview.isSelected();
        updatePreview();
    }//GEN-LAST:event_jcbPreviewActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    protected javax.swing.JCheckBox jcbPreview;
    protected browser.filebrowsers.JListFileFilter jlFileFilter;
    protected javax.swing.JLabel jlFiltering;
    protected javax.swing.JPanel jpCenter;
    protected javax.swing.JPanel jpFileBrowser;
    browser.filebrowsers.JPanelImageInfo jpImageInfo;
    protected javax.swing.JPanel jpPreview;
    private javax.swing.JScrollPane jspFilesFiltering;
    // End of variables declaration//GEN-END:variables
}
