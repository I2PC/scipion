/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.Cache;
import browser.ICONS_MANAGER;
import browser.LABELS;
import browser.filebrowsers.JPanelXmippBrowser;
import browser.imageitems.listitems.AbstractImageItem;
import browser.imageitems.listitems.FileItem;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Toolbar;
import ij.measure.Calibration;
import java.awt.Image;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JLabel;

/**
 *
 * @author Juanjo Vega
 */
abstract public class JPanelXmippFilter extends JPanelXmippBrowser {

    JLabel jlFilter;
    Cache<String, ImagePlus> cache = new Cache<String, ImagePlus>();
    final static int W = 512, H = 512;

    public JPanelXmippFilter(String folder) {
        this(folder, "");
    }

    public JPanelXmippFilter(String folder, String expression) {
        super(folder, expression);

        setPreviewSize(256 + 10, 256 + 10);

        jpCenter.remove(jlFiltering);   // Hide filter stuff.
        jpFileBrowser.remove(searchBox);
        jpPreview.remove(jcbPreview);

        jlFilter = new JLabel();
        jlFilter.setBorder(BorderFactory.createTitledBorder(LABELS.TITLE_PREVIEW_FILTER));
        jlFilter.setPreferredSize(jpPreview.getPreferredSize());

        jpPreview.add(jlFilter, java.awt.BorderLayout.NORTH);
    }

    @Override
    protected void setPreview(AbstractImageItem imageItem) {
        System.out.println("super: Set preview");
        super.setPreview(imageItem);

        System.out.println("super: Set FILTERED preview");
        Image filteredPreview = null;

        try {
            filteredPreview = getFilteredPreview(imageItem).getImage();
        } catch (Exception e) {
            filteredPreview = ICONS_MANAGER.MISSING_ITEM.getImage();
        }

        System.out.println("Fin de set filter preview");
        jlFilter.setIcon(new ImageIcon(filteredPreview));
    }

    @Override
    protected void clearPreview() {
        super.clearPreview();
        jlFilter.setIcon(null);
    }

    abstract ImagePlus getFilteredPreview(AbstractImageItem item) throws Exception;

    @Override
    protected void openFileAsDefault(FileItem item) {
        try {
            if (item instanceof AbstractImageItem) {
                ImagePlus preview = getFilteredPreview((AbstractImageItem) item);

                ImagePlus imp = new ImagePlus(preview.getTitle(),
                        preview.getProcessor().resize(W, H));

                Calibration c = new Calibration();
                c.pixelWidth = c.pixelHeight = 1. / imp.getWidth();
                imp.setCalibration(c);

                // Check if ImageJ is opened.
                if (IJ.getInstance() == null) {
                    ImageJ ij = new ImageJ();
                }
                // ...to set line tool beforehand.
                IJ.setTool(Toolbar.LINE);

                ImagesWindowFactory.captureFrame(imp);
            }
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
        }
    }

    @Override
    protected void openSelectedFile() {
        //if (jlFileFilter.getSelectedIndex() > 0) {  // Avoid parent...
        FileItem item = (FileItem) jlFileFilter.getSelectedValue();

        if (!item.isDirectory()) {  // ...and directories.
            openFileAsDefault(item);
        }
        //}
    }

    @Override
    public void refreshCurrentDirectory() {
        super.refreshCurrentDirectory();
        cache.clear();
    }
}
