
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
        return convertToImagej(image, path, true);
    }

    public static ImagePlus convertToImagej(ImageDouble image, String path, boolean useLogarithm) {
        if (image.isPSD()) {
            image.convertPSD(useLogarithm);
        }

        int w = image.getXsize();
        int h = image.getYsize();
        int d = image.getZsize();
        long n = image.getNsize();

        return convertToImagej(image.getData(), w, h, d, n, path);
    }

    private static ImagePlus convertToImagej(double array[], int w, int h, int d, long n, String title) {
        int sliceSize = w * h;
        int imageSize = sliceSize * d;
        double out[] = new double[sliceSize];
        ImageStack is = new ImageStack(w, h);

        for (int i = 0; i < n; i++) {
            int offset = i * imageSize;
            for (int j = 0; j < d; j++) {
                System.arraycopy(array, offset + j * sliceSize, out, 0, sliceSize);

                FloatProcessor processor = new FloatProcessor(w, h, out);
                is.addSlice(String.valueOf(i), processor);
            }
        }

        return new ImagePlus(title, is);
    }
}
