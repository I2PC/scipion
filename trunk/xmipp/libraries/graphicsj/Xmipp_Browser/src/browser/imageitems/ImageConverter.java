package browser.imageitems;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import xmipp.ImageDouble;
import xmipp.MultidimArrayd;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class ImageConverter {

    public static ImagePlus xmipp2Imagej(ImageDouble img, String filename) {
        if (img.getData() == null) {
            return new ImagePlus();
        }

        int w = img.getData().getXdim();
        int h = img.getData().getYdim();

        // Creates image
        FloatProcessor ip = new FloatProcessor(w, h);

        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                ip.setf(i, j, (float) img.getPixel(j, i));
            }
        }

        ImagePlus imagePlus = new ImagePlus(filename, ip);

        return imagePlus;
    }

    public static MultidimArrayd imagej2xmipp(ImagePlus image) {
        int w = image.getWidth();
        int h = image.getHeight();
        int d = image.getStackSize();

        MultidimArrayd matrix = new MultidimArrayd(d, h, w);

        for (int k = 0; k < d; k++) {
            ImageProcessor processor = image.getStack().getProcessor(k + 1);
            for (int i = 0; i < h; i++) {
                for (int j = 0; j < w; j++) {
                    matrix.setVoxel(k, i, j, processor.getPixelValue(j, i));
                }
            }
        }

        matrix.setXmippOrigin();

        return matrix;
    }

    public static ImagePlus xmipp2imagej(MultidimArrayd matrix) {
        int w = matrix.getXdim();
        int h = matrix.getYdim();
        int d = matrix.getZdim();

        float image[][] = new float[d][w * h];

        for (int k = 0; k < d; k++) {
            int pointer = 0;
            for (int j = 0; j < h; j++) {
                for (int i = 0; i < w; i++, pointer++) {
                    image[k][pointer] = (float) matrix.getVoxel(k, i, j);
                }
            }
        }

        // Builds the stack.
        ImageStack outputStack = new ImageStack(w, h, d);
        for (int slice = 0; slice < image.length; slice++) {
            outputStack.setPixels(image[slice], slice + 1);
        }

        return new ImagePlus("", outputStack);
    }
}
