/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.filebrowsers;

import browser.LABELS;
import browser.SpringUtilities;
import browser.imageitems.listitems.AbstractImageItem;
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
public class JPanelXmippFileListPSD extends JPanelXmippBrowser {

    protected JLabel jlDownsampling = new JLabel(LABELS.LABEL_DOWNSAMPLING);
    protected JSpinner jsDownsampling = new JSpinner();

    public JPanelXmippFileListPSD(String folder) {
        this(folder, "");
    }

    public JPanelXmippFileListPSD(String folder, String expression) {
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

        try {
            ImageDouble img = new ImageDouble();
            img.readHeader(filename);

            double downsampling = (Integer) jsDownsampling.getValue();
            double pixels[] = ImageDouble.fastEstimateEnhancedPSD(filename, downsampling);

            FloatProcessor ip = new FloatProcessor(img.getXsize(), img.getYsize(), pixels);
            ImagePlus imp = new ImagePlus(filename, ip);

            // @TODO COSS: Scale ? (*)
            //imp.getProcessor().scale(ICONS_MANAGER.DEFAULT_PREVIEW_WIDTH, ICONS_MANAGER.DEFAULT_PREVIEW_HEIGHT);
            imp.show();
            // return imp;
        } catch (Exception ex) {
            IJ.error(filename + " >>> " + ex.getMessage());
        }

        return null;
    }
}
