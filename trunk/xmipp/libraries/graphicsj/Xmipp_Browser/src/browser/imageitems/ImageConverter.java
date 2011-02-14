package browser.imageitems;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import xmipp.ImageDouble;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class ImageConverter {

    public static ImagePlus convertToImagej(ImageDouble image, String path) {
        if (image.isPSD()) {
            image.convertPSD();
        }

        int w = image.getWidth();
        int h = image.getHeight();
        int d = image.getDepth();
        long n = image.getNimages();

        return convertToImagej(image.getData(), w, h, d, path);
    }

    private static ImagePlus convertToImagej(double array[], int w, int h, int d, String title) {
        double out[][] = new double[d][w * h];

        // @TODO ArrayCopy
        for (int k = 0; k < d; k++) {
            for (int i = 0; i < h; i++) {
                for (int j = 0; j < w; j++) {
                    out[k][j * w + i] = array[j * w + i];
                }
            }
        }

        ImageStack is = new ImageStack(w, h);
        for (int i = 0; i < out.length; i++) {
            FloatProcessor processor = new FloatProcessor(w, h, out[i]);
            is.addSlice(String.valueOf(i), processor);
        }

        return new ImagePlus(title, is);
    }
}
