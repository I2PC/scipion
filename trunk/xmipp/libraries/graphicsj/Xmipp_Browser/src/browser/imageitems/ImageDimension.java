/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems;

import ij.IJ;
import xmipp.ImageGeneric;

/**
 *
 * @author Juanjo Vega
 */
public class ImageDimension {

    private int width, height, depth;
    private long nimages;

    public ImageDimension() {
    }

    public ImageDimension(int width, int height, int depth, long nimages) {
        this.width = width;
        this.height = height;
        this.depth = depth;
        this.nimages = nimages;
    }

    public ImageDimension(ImageGeneric image) {
        try {
            width = image.getXDim();
            height = image.getYDim();
            depth = image.getZDim();
            nimages = image.getNDim();
        } catch (Exception ex) {
            IJ.error("Retrieving image dimensions: "+ image);
        }
//        System.out.println(nimages + " images.");
    }

    public int getWidth() {
        return width;
    }

    public int getHeight() {
        return height;
    }

    public int getDepth() {
        return depth;
    }

    public long getNimages() {
        return nimages;
    }

    public void setWidth(int width) {
        this.width = width;
    }

    public void setHeight(int height) {
        this.height = height;
    }

    public void setDepth(int depth) {
        this.depth = depth;
    }

    public void setNImages(long nimages) {
        this.nimages = nimages;
    }
}
