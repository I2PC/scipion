
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

//        System.out.println("(" + w + ", " + h + ", " + d + ", " + n + ")");

        return convertToImagej(image.getData(), w, h, d, path);
    }

    private static ImagePlus convertToImagej(double array[], int w, int h, int d, String title) {
        int sliceSize = w * h;
        double out[] = new double[sliceSize];
        ImageStack is = new ImageStack(w, h);

        for (int i = 0; i < d; i++) {
            System.arraycopy(array, i * sliceSize, out, 0, sliceSize);

            FloatProcessor processor = new FloatProcessor(w, h, out);
            is.addSlice(String.valueOf(i), processor);
        }

        return new ImagePlus(title, is);
    }
}
