/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.LABELS;
import browser.imageitems.ImageConverter;
import browser.imageitems.listitems.AbstractImageItem;
import ij.IJ;
import ij.ImagePlus;
import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public class JPanelXmippBadPixelsFilter extends JPanelXmippGaussianFilter {

    public JPanelXmippBadPixelsFilter(String metadata, double bad_pixels_factor) {
        super(metadata, bad_pixels_factor);

        jlW1.setText(LABELS.LABEL_BAD_PIXELS_FACTOR);
    }

    public double getFactor() {
        return getW1();
    }

    @Override
    ImagePlus getFilteredPreview(AbstractImageItem item) {
        ImagePlus imp = null;

        try {
            double factor = (Double) jsW1.getValue();
            ImageDouble image = new ImageDouble();
            image.badPixelsFilter(item.getAbsoluteFileName(),
                    factor, previewWidth, previewHeight);

            imp = ImageConverter.convertToImageJ(image, item.getFileName());
        } catch (Exception e) {
            IJ.error(e.getMessage());
        }

        return imp;
    }
}
