/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.filebrowsers;

import browser.Cache;
import browser.ICONS_MANAGER;
import browser.LABELS;
import browser.SpringUtilities;
import browser.imageitems.listitems.AbstractImageItem;
import browser.imageitems.listitems.FileItem;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import java.awt.BorderLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpringLayout;
import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public class JPanelXmippFileListCTF extends JPanelXmippBrowser {

    protected JLabel jlDownsampling = new JLabel(LABELS.LABEL_DOWNSAMPLING);
    protected JSpinner jsDownsampling = new JSpinner();
    Cache<String, ImagePlus> cache = new Cache<String, ImagePlus>();

    public JPanelXmippFileListCTF(String folder) {
        this(folder, "");
    }

    public JPanelXmippFileListCTF(String folder, String expression) {
        super(folder, expression);

        jsDownsampling.setValue(2);

        JPanel jpDownsampling = new JPanel(new SpringLayout());

        jpDownsampling.add(jlDownsampling);
        jpDownsampling.add(jsDownsampling);

        SpringUtilities.makeCompactGrid(jpDownsampling, 1, 2, 3, 3, 3, 3);

        jpCenter.add(jpDownsampling, BorderLayout.NORTH);
    }

    @Override
    protected ImagePlus getPreview(AbstractImageItem item) {
        String filename = item.getAbsoluteFileName();

        ImagePlus imp = (ImagePlus) cache.get(filename);

        if (imp == null) {
            try {
                double downsampling = (Integer) jsDownsampling.getValue();
                double pixels[] = ImageDouble.fastEstimateEnhancedPSD(
                        filename,
                        downsampling,
                        ICONS_MANAGER.DEFAULT_PREVIEW_WIDTH,
                        ICONS_MANAGER.DEFAULT_PREVIEW_HEIGHT);

                FloatProcessor ip = new FloatProcessor(ICONS_MANAGER.DEFAULT_PREVIEW_WIDTH,
                        ICONS_MANAGER.DEFAULT_PREVIEW_HEIGHT, pixels);
                imp = new ImagePlus(filename, ip);

                cache.put(filename, imp);
            } catch (Exception ex) {
                IJ.error(ex.getMessage());
            }
        }

        return imp;
    }

    @Override
    protected void openFileAsDefault(FileItem item) {
        String filename = item.getAbsoluteFileName();
        double downsampling = (Integer) jsDownsampling.getValue();

        try {
            ImageDouble img = new ImageDouble();
            img.read(filename);

            double pixels[] = ImageDouble.fastEstimateEnhancedPSD(filename, downsampling);

            FloatProcessor ip = new FloatProcessor(img.getXsize(), img.getYsize(), pixels);
            ImagePlus imp = new ImagePlus(filename, ip);

            ImagesWindowFactory.captureFrame(imp);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
        }
    }

    @Override
    public void refreshCurrentDirectory() {
        super.refreshCurrentDirectory();
        cache.clear();
    }
}
