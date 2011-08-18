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
import ij.gui.Toolbar;
import ij.measure.Calibration;
import ij.process.FloatProcessor;
import java.awt.BorderLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SpringLayout;
import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public class JPanelXmippFileListPSD extends JPanelXmippBrowser {

    protected JLabel jlDownsampling = new JLabel(LABELS.LABEL_DOWNSAMPLING);
    protected JSpinner jsDownsampling = new JSpinner();
    Cache<String, ImagePlus> cache = new Cache<String, ImagePlus>();
//    private final static int DELAY_TO_UPDATE = 500;
//    private Timer updateTimer = new Timer(true);    // Timer for preview update.
//    private PreviewUpdater previewUpdaterTask;    // Associated task for preview update.

    public JPanelXmippFileListPSD(String folder) {
        this(folder, "");
    }

    public JPanelXmippFileListPSD(String folder, String expression) {
        super(folder, expression);

        SpinnerNumberModel model = new SpinnerNumberModel();
        model.setStepSize(new Double(0.1));
        model.setValue(new Double(2.0));
        jsDownsampling = new JSpinner(model);
        jsDownsampling.setEditor(new JSpinner.NumberEditor(jsDownsampling, "#.###"));
        ((JSpinner.DefaultEditor) jsDownsampling.getEditor()).getTextField().setColumns(3);
        jsDownsampling.addChangeListener(new javax.swing.event.ChangeListener() {

            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                cache.clear();
                updatePreview();
                //startUpdater();
            }
        });


        JPanel jpDownsampling = new JPanel(new SpringLayout());

        jpDownsampling.add(jlDownsampling);
        jpDownsampling.add(jsDownsampling);

        jpCenter.remove(jlFiltering);   // Hide filter stuff.
        jpFileBrowser.remove(searchBox);
        jpPreview.remove(jcbPreview);

        SpringUtilities.makeCompactGrid(jpDownsampling, 1, 2, 3, 3, 3, 3);

        jpCenter.add(jpDownsampling, BorderLayout.NORTH);
    }

    public double getDownsampling() {
        return (Double) jsDownsampling.getValue();
    }
//
//    private void startUpdater() {
//        if (previewUpdaterTask != null) {
//            previewUpdaterTask.cancel();
//        }
//
//        previewUpdaterTask = new PreviewUpdater();
//        updateTimer.schedule(previewUpdaterTask, DELAY_TO_UPDATE);
//    }

    @Override
    protected ImagePlus getPreview(AbstractImageItem item) {
        String filename = item.getAbsoluteFileName();

        ImagePlus imp = (ImagePlus) cache.get(filename);

        if (imp == null) {
            try {
                double downsampling = (Double) jsDownsampling.getValue();
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
    protected void openSelectedFile() {
        if (jlFileFilter.getSelectedIndex() > 0) {  // Avoid parent...
            FileItem item = (FileItem) jlFileFilter.getSelectedValue();

            if (!item.isDirectory()) {  // ...and directories.
                openFileAsDefault(item);
            }
        }
    }

    @Override
    protected void openFileAsDefault(FileItem item) {
        String filename = item.getAbsoluteFileName();
        double downsampling = (Double) jsDownsampling.getValue();

        try {
            ImageDouble img = new ImageDouble();
            img.readHeader(filename);

            double pixels[] = ImageDouble.fastEstimateEnhancedPSD(
                    filename, downsampling, img.getXsize(), img.getYsize());

            FloatProcessor ip = new FloatProcessor(img.getXsize(), img.getYsize(), pixels);
            ImagePlus imp = new ImagePlus(filename, ip);

            Calibration c = new Calibration();
            c.pixelWidth = c.pixelHeight = 1. / ip.getWidth();
            imp.setCalibration(c);

            // Sets line tool beforehand.
            IJ.setTool(Toolbar.LINE);

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
//
//    class PreviewUpdater extends TimerTask {
//
//        @Override
//        public void run() {
//            cache.clear();
//            updatePreview();
//        }
//    }
}
