/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.filebrowsers;

import browser.Cache;
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
public class JPanelXmippFileListCTF extends JPanelXmippBrowser {

    protected JLabel jlDownsampling = new JLabel(LABELS.LABEL_DOWNSAMPLING);
    protected JSpinner jsDownsampling = new JSpinner();
    Cache<String, ImagePlus> cache = new Cache<String, ImagePlus>();
    final static int W = 512, H = 512;
//    private final static int DELAY_TO_UPDATE = 500;
//    private Timer updateTimer = new Timer(true);    // Timer for preview update.
//    private PreviewUpdater previewUpdaterTask;    // Associated task for preview update.

    public JPanelXmippFileListCTF(String folder) {
        this(folder, "");
    }

    public JPanelXmippFileListCTF(String folder, String expression) {
        this(folder, expression, 1.0);
    }

    public JPanelXmippFileListCTF(String folder, String expression, double downsampling) {
        super(folder, expression);

        setPreviewSize(256 + 10, 256 + 10);

        SpinnerNumberModel model = new SpinnerNumberModel();
        model.setStepSize(new Double(0.1));
        model.setValue(new Double(downsampling));
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

    public void setDownsamplingEnabled(boolean enabled) {
        jsDownsampling.setEnabled(enabled);
    }

    public double getDownsampling() {
        return (Double) jsDownsampling.getValue();
    }

    @Override
    protected ImagePlus getPreview(AbstractImageItem item) {
        String filename = item.getAbsoluteFileName();

        ImagePlus imp = (ImagePlus) cache.get(filename);

        if (imp == null) {
            try {
                double downsampling = (Double) jsDownsampling.getValue();
                double pixels[] = ImageDouble.fastEstimateEnhancedPSD(
                        filename, downsampling, previewWidth, previewHeight);

                FloatProcessor ip = new FloatProcessor(
                        previewWidth, previewHeight, pixels);
                imp = new ImagePlus(item.getFileName(), ip);

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
        if (item instanceof AbstractImageItem) {
            ImagePlus preview = getPreview((AbstractImageItem) item);

            ImagePlus imp = new ImagePlus(preview.getTitle(),
                    preview.getProcessor().resize(W, H));

            Calibration c = new Calibration();
            c.pixelWidth = c.pixelHeight = 1. / imp.getWidth();
            imp.setCalibration(c);

            // Sets line tool beforehand.
            IJ.setTool(Toolbar.LINE);

            ImagesWindowFactory.captureFrame(imp);
        }
    }
//    protected void openFileAsDefault(FileItem item) {
//        String filename = item.getAbsoluteFileName();
//        double downsampling = (Double) jsDownsampling.getValue();
//
//        try {
//            ImageDouble img = new ImageDouble();
//            img.readHeader(filename);
//
//            double pixels[] = ImageDouble.fastEstimateEnhancedPSD(
//                    filename, downsampling, img.getXsize(), img.getYsize());
//
//            FloatProcessor ip = new FloatProcessor(img.getXsize(), img.getYsize(), pixels);
//            ImagePlus imp = new ImagePlus(filename, ip);
//
//            Calibration c = new Calibration();
//            c.pixelWidth = c.pixelHeight = 1. / ip.getWidth();
//            imp.setCalibration(c);
//
//            // Sets line tool beforehand.
//            IJ.setTool(Toolbar.LINE);
//
//            ImagesWindowFactory.captureFrame(imp);
//        } catch (Exception ex) {
//            IJ.error(ex.getMessage());
//        }
//    }

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
