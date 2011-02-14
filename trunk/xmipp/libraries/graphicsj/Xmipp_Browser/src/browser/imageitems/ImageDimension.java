/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems;

import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public class ImageDimension {

    public int width, height, depth;
    public long nimages;

    public ImageDimension() {
    }

    public ImageDimension(int width, int height, int depth, long nimages) {
        this.width = width;
        this.height = height;
        this.depth = depth;
        this.nimages = nimages;
    }

    public ImageDimension(ImageDouble image) {
        width = image.getWidth();
        height = image.getHeight();
        depth = image.getDepth();
        nimages = image.getNimages();

//        System.out.println(nimages + " images.");
    }
}
