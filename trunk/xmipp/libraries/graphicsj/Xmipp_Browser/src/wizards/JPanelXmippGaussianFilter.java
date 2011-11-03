/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.LABELS;
import browser.SpringUtilities;
import browser.imageitems.ImageConverter;
import browser.imageitems.listitems.AbstractImageItem;
import ij.IJ;
import ij.ImagePlus;
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
public class JPanelXmippGaussianFilter extends JPanelXmippFilterMetadata {

    JPanel jpSpinners;
    JLabel jlW1 = new JLabel(LABELS.LABEL_W1);
    JSpinner jsW1 = new JSpinner();

    public JPanelXmippGaussianFilter(String metadata) {
        this(metadata, 1.0);
    }

    public JPanelXmippGaussianFilter(String metadata, double w1) {
        super(metadata);

        SpinnerNumberModel model = new SpinnerNumberModel();
        model.setStepSize(new Double(0.1));
        model.setValue(new Double(w1));
        jsW1 = new JSpinner(model);
        jsW1.setEditor(new JSpinner.NumberEditor(jsW1, "#.###"));
        ((JSpinner.DefaultEditor) jsW1.getEditor()).getTextField().setColumns(3);
        jsW1.addChangeListener(new javax.swing.event.ChangeListener() {

            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                cache.clear();
                updatePreview();
                //startUpdater();
            }
        });

        jpSpinners = new JPanel(new SpringLayout());
        jpSpinners.add(jlW1);
        jpSpinners.add(jsW1);

        jpCenter.remove(jlFiltering);   // Hide filter stuff.
        jpFileBrowser.remove(searchBox);
        jpPreview.remove(jcbPreview);

        SpringUtilities.makeCompactGrid(jpSpinners, 1, 2, 3, 3, 3, 3);

        jpCenter.add(jpSpinners, BorderLayout.NORTH);
    }

    public double getW1() {
        return (Double) jsW1.getValue();
    }

    @Override
    ImagePlus getFilteredPreview(AbstractImageItem item) {
        ImagePlus imp = null;

        try {
            double w1 = (Double) jsW1.getValue();
            ImageDouble image = new ImageDouble();
            image.gaussianFilter(item.getAbsoluteFileName(),
                    w1, previewWidth, previewHeight);

            imp = ImageConverter.convertToImageJ(image, item.getFileName());
        } catch (Exception e) {
            IJ.error(e.getMessage());
        }

        return imp;
    }
}
