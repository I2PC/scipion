
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class ImageConverter {

    public static ImagePlus xmipp2imagej(Projection projection, int w, int h) {
        float image[] = new float[w * h];
        int i0 = h / 2;
        int j0 = w / 2;

        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                image[(h - 1 - i) * w + w - 1 - j] = (float) projection.getPixel(i - i0, j - j0);
            }
        }

        // Creates image.
        FloatProcessor ip = new FloatProcessor(w, h);
        ImagePlus imagePlus = new ImagePlus("Projection", ip);
        imagePlus.getProcessor().setPixels(image);

        return imagePlus;
    }

    /**
     * Converts xmipp into imageplus.
     * @param matrix
     * @param w
     * @param h
     * @param d
     * @return
     */
    public static ImagePlus xmipp2imagej_(MultidimArrayd matrix, int w, int h, int d) {
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

    /**
     * Converts imageplus into xmipp.
     * @param image
     * @return
     */
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

    /**
     * Converts imageplus into xmipp projection (2D image).
     * @param image
     * @return
     */
    public static Projection imagej2xmippProjection(ImagePlus image) {
        int w = image.getWidth();
        int h = image.getHeight();


        Projection projection = new Projection();
        projection.reset(h, w);

        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                projection.setPixel(i, j, image.getProcessor().getPixelValue(j, i));
            }
        }

        return projection;
    }
}
