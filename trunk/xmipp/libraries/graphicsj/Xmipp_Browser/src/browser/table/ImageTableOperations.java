/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table;

import browser.imageitems.TableImageItem;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import java.util.Vector;

/**
 *
 * @author Juanjo Vega
 */
public class ImageTableOperations {

    public static ImagePlus mean(ImagePlus images[]) {
        int w = images[0].getWidth();
        int h = images[0].getHeight();
        float mean[][] = new float[w][h];

        // Traverses image size.
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                float pixel = 0;
                for (int k = 0; k < images.length; k++) {   // Deepness
                    pixel += images[k].getProcessor().getPixelValue(i, j);
                }

                mean[i][j] = pixel / images.length;    // Sets mean value.
            }
        }

        FloatProcessor processor = new FloatProcessor(mean);
        ImagePlus meanIp = new ImagePlus();
        meanIp.setProcessor("", processor);

        return meanIp;
    }

    public static ImagePlus std_deviation(ImagePlus images[]) {
        int w = images[0].getWidth();
        int h = images[0].getHeight();

        ImagePlus mean = mean(images);

        float std[] = new float[w * h];

        // Traverses image's size.
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                for (int k = 0; k < images.length; k++) {   // Deepness
                    float pixel = images[k].getProcessor().getPixelValue(i, j);
                    float avg = mean.getProcessor().getPixelValue(i, j);

//                    sum2 += (pixel - avg) * (pixel - avg);
                    std[j * w + i] += (pixel - avg) * (pixel - avg);
                }
                std[j * w + i] /= images.length;
            }
        }

        FloatProcessor processor = new FloatProcessor(w, h);
        processor.setPixels(std);
        ImagePlus std_ip = new ImagePlus();
        std_ip.setProcessor("", processor);

        return std_ip;
    }

    public static double[] getMinAndMax(Vector<TableImageItem> items) {
        // Calculates min and max.
        double min = items.elementAt(0).getPreview().getProcessor().getMin();
        double max = items.elementAt(0).getPreview().getProcessor().getMax();

        for (int i = 1; i < items.size(); i++) {
            double current_min = items.elementAt(i).getPreview().getProcessor().getMin();
            double current_max = items.elementAt(i).getPreview().getProcessor().getMax();

            if (current_min < min) {
                min = current_min;
            }

            if (current_max > max) {
                max = current_max;
            }
        }

        return new double[]{min, max};
    }

    public static void normalize(Vector<TableImageItem> items) {
        // Calculates min and max.
        double min = items.elementAt(0).getPreview().getProcessor().getMin();
        double max = items.elementAt(0).getPreview().getProcessor().getMax();

        for (int i = 1; i < items.size(); i++) {
            double current_min = items.elementAt(i).getPreview().getProcessor().getMin();
            double current_max = items.elementAt(i).getPreview().getProcessor().getMax();

            if (current_min < min) {
                min = current_min;
            }

            if (current_max > max) {
                max = current_max;
            }
        }

        // Calculates factor.
        double factor = 255 / (max - min);

        System.out.println(" >>> Max: " + max + " / Min: " + min + " / Factor: " + factor);

        // Apllies factor to normalized images array.
        for (int i = 0; i < items.size(); i++) {
            items.elementAt(i).setNormalized(factor, min);

            // Show normalized.
/*            ImagePlus imp = items.elementAt(i).getImagePlus();
            normalize(imp, factor, min);
            ImagesWindowFactory.openImageWindow(imp);
            ImagesWindowFactory.openImageWindow(items.elementAt(i).getFile());*/
        }
    }

    /**
     * Normalizes image by: factor * (I(x,y) - min(I))
     * factor = #levels / (max(I) - min(I))
     * @param image
     * @param factor
     * @param min
     */
    /*public static void normalize(ImagePlus image, double factor, double min) {
    double current_min = image.getProcessor().getMin();
    double current_max = image.getProcessor().getMax();

    double new_min = (current_min - min) * factor;
    double new_max = (current_max - min) * factor;

    image.getProcessor().setMinAndMax(new_min, new_max);
    image.updateImage();//updateAndDraw();
    }*/
    /*
    public static void normalize_(ImagePlus image, double factor, double min) {
    float pixels[] = (float[]) image.getProcessor().getPixels();
    float new_pixels[] = new float[pixels.length];

    int w = image.getWidth();
    int h = image.getHeight();

    // Traverses image size.
    for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
    //                System.out.print(" -> " + pixels[j * w + i]);

    new_pixels[j * w + i] = (float) (factor * ((pixels[j * w + i]) - min));

    //                System.out.println(" >> " + pixels[j * w + i]);
    }
    }

    image.getProcessor().setPixels(new_pixels);
    image.updateImage();//updateAndDraw();
    }*/
    public static ImagePlus clone(ImagePlus image) {
        float pixels[][] = new float[image.getWidth()][image.getHeight()];

        // Traverses image size.
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixels[i][j] = image.getProcessor().getPixelValue(i, j);
            }
        }

        FloatProcessor processor = new FloatProcessor(pixels);
        ImagePlus normalized = new ImagePlus();
        normalized.setProcessor("", processor);

        return normalized;
    }
}
