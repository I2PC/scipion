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

    /*
    public static void main(String args[]) {
    int w, h;
    w = h = 5;
    double input[] = new double[w * h];

    for (int i = 0; i < w; i++) {
    for (int j = 0; j < h; j++) {
    input[j * w + i] = Double.parseDouble(j + "." + i);
    }
    }

    print(input, w, h);
    double out[] = traverse(input, w, h);
    System.out.println(" -- -- -- -- -- -- -- -- ");
    print(out, w, h);
    }

    public static double[] traverse(double in[], int w, int h) {
    double out[] = new double[w * h];

    for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
    out[j * w + i] = in[j * w + i];
    }
    }

    return out;
    }

    public static void print(double array[], int w, int h) {
    for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
    System.out.print(" " + array[j * w + i]);
    }
    System.out.println();
    }
    }*/
}
