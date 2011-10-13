/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.LABELS;
import browser.SpringUtilities;
import browser.imageitems.ImageConverter;
import browser.imageitems.listitems.AbstractImageItem;
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
public class JPanelXmippFileListCTF extends JPanelXmippFilter {

    protected JLabel jlDownsampling = new JLabel(LABELS.LABEL_DOWNSAMPLING);
    protected JSpinner jsDownsampling = new JSpinner();

    public JPanelXmippFileListCTF(String folder) {
        this(folder, "");
    }

    public JPanelXmippFileListCTF(String folder, String expression) {
        this(folder, expression, 1.0);
    }

    public JPanelXmippFileListCTF(String folder, String expression, double downsampling) {
        super(folder, expression);

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
    ImagePlus getFilteredPreview(AbstractImageItem item) throws Exception {
        double downsampling = (Double) jsDownsampling.getValue();
        ImageDouble image = new ImageDouble();
        image.fastEstimateEnhancedPSD(item.getAbsoluteFileName(),
                downsampling, previewWidth, previewHeight);

        return ImageConverter.convertToImageJ(image, item.getFileName());
    }
}
