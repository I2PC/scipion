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
import ij.process.FloatProcessor;
import javax.swing.JLabel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public class JPanelXmippBandPassFilter extends JPanelXmippGaussianFilter {

    protected JLabel jlW2 = new JLabel(LABELS.LABEL_W2);
    protected JLabel jlRaisedW = new JLabel(LABELS.LABEL_RAISED_W);
    protected JSpinner jsW2 = new JSpinner();
    protected JSpinner jsRaisedW = new JSpinner();

    public JPanelXmippBandPassFilter(String metadata) {
        this(metadata, 1.0, 1.0, 1.0);
    }

    public JPanelXmippBandPassFilter(String metadata, double w1, double w2, double raised_w) {
        super(metadata, w1);

        SpinnerNumberModel model = new SpinnerNumberModel();
        model.setStepSize(new Double(0.1));
        model.setValue(new Double(w2));
        jsW2 = new JSpinner(model);
        jsW2.setEditor(new JSpinner.NumberEditor(jsW2, "#.###"));
        ((JSpinner.DefaultEditor) jsW2.getEditor()).getTextField().setColumns(3);
        jsW2.addChangeListener(new javax.swing.event.ChangeListener() {

            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                cache.clear();
                updatePreview();
                //startUpdater();
            }
        });

        model = new SpinnerNumberModel();
        model.setStepSize(new Double(0.1));
        model.setValue(new Double(raised_w));
        jsRaisedW = new JSpinner(model);
        jsRaisedW.setEditor(new JSpinner.NumberEditor(jsRaisedW, "#.###"));
        ((JSpinner.DefaultEditor) jsW2.getEditor()).getTextField().setColumns(3);
        jsRaisedW.addChangeListener(new javax.swing.event.ChangeListener() {

            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                cache.clear();
                updatePreview();
                //startUpdater();
            }
        });

        jpSpinners.add(jlW2);
        jpSpinners.add(jsW2);
        jpSpinners.add(jlRaisedW);
        jpSpinners.add(jsRaisedW);

        SpringUtilities.makeCompactGrid(jpSpinners, 3, 2, 3, 3, 3, 3);
    }

    public double getW2() {
        return (Double) jsW2.getValue();
    }

    public double getRaisedW() {
        return (Double) jsRaisedW.getValue();
    }

    @Override
    ImagePlus getFilteredPreview(AbstractImageItem item) {
        ImagePlus imp = null;

        try {
            double w1 = (Double) jsW1.getValue();
            double w2 = (Double) jsW2.getValue();
            double raised_w = (Double) jsRaisedW.getValue();
            ImageDouble image = new ImageDouble();
            image.bandPassFilter(item.getAbsoluteFileName(),
                    w1, w2, raised_w,
                    previewWidth, previewHeight);

            imp = ImageConverter.convertToImageJ(image, item.getFileName());
        } catch (Exception e) {
            IJ.error(e.getMessage());
        }

        return imp;
    }
}
