package xmipp.io.writers;

import ij.ImagePlus;
import xmippij.XmippImageConverter;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class ImageWriter extends Writer {

    @Override
    public void write(ImagePlus imp, String path) throws Exception {
        XmippImageConverter.convertToXmipp(imp).write(path);
    }

    @Override
    protected String getSaveDialogTitle() {
        return "Save as xmipp image...";
    }
}
