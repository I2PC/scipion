/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.ij.plugins.maskstoolbar;

import ij.ImagePlus;

/**
 *
 * @author Juanjo Vega
 */
public class Operations {

    public static void binarize(ImagePlus ip) {
        byte pixels[] = (byte[]) ip.getProcessor().getPixels();

        for (int i = 0; i < pixels.length; i++) {
            if (pixels[i] > 0) {
                pixels[i] = 1;
            }
        }

        ip.getProcessor().setPixels(pixels);
        //ip.updateImage();
    }
}
