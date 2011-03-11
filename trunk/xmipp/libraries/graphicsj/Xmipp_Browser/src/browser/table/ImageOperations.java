/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table;

import browser.imageitems.TableImageItem;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.util.Vector;

/**
 *
 * @author Juanjo Vega
 */
public class ImageOperations {

    public static ImagePlus average(Vector<TableImageItem> items) {
        ImagePlus current = items.firstElement().getImagePlus();
        int w = current.getWidth();
        int h = current.getHeight();
        double average[] = new double[w * h];

        // For all images...
        for (int k = 0; k < items.size(); k++) {
            current = items.elementAt(k).getImagePlus();

            // Adds current image to sum.
            for (int j = 0; j < h; j++) {
                for (int i = 0; i < w; i++) {
                    average[j * w + i] += current.getProcessor().getPixelValue(i, j);
                }
            }

            current.close();
            current = items.elementAt(k).getImagePlus();
        }

        // Calculates average...
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                average[j * w + i] /= items.size(); // ...by dividing with the #images
            }
        }

        FloatProcessor processor = new FloatProcessor(w, h, average);
        ImagePlus averageIp = new ImagePlus();
        averageIp.setProcessor("", processor);

        return averageIp;
    }

    public static ImagePlus std_deviation(Vector<TableImageItem> items) {
        ImagePlus current = items.firstElement().getImagePlus();
        int w = current.getWidth();
        int h = current.getHeight();

        ImagePlus average = average(items);
        double std[] = new double[average.getWidth() * average.getHeight()];

        // Traverses images.
        for (int k = 0; k < items.size(); k++) {
            current = items.elementAt(k).getImagePlus();

            for (int j = 0; j < h; j++) {
                for (int i = 0; i < w; i++) {
                    float pixel = current.getProcessor().getPixelValue(i, j);
                    float avg = average.getProcessor().getPixelValue(i, j);

                    std[j * w + i] += (pixel - avg) * (pixel - avg);
                }
            }
            current.close();
        }

        for (int i = 0; i < h; i++) {
            for (int j = 0; j < h; j++) {
                //std[j * w + i] /= items.size();
                std[j * w + i] = Math.sqrt(std[j * w + i] / items.size());
            }
        }

        FloatProcessor processor = new FloatProcessor(w, h, std);
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

    public static double mean(ImagePlus imp) {
        double mean = 0.0;

        FloatProcessor fp = new FloatProcessor(imp.getWidth(), imp.getHeight());
        ImageProcessor ip = imp.getProcessor();
        ip.toFloat(0, fp);

        float pixels[] = (float[]) fp.getPixels();
        for (int i = 0; i < pixels.length; i++) {
            mean += pixels[i];
        }

        mean /= pixels.length;

        return mean;
    }

    public static double std_dev(ImagePlus imp) {
        double mean = mean(imp);
        double std_dev = 0.0;

        FloatProcessor fp = new FloatProcessor(imp.getWidth(), imp.getHeight());
        ImageProcessor ip = imp.getProcessor();
        ip.toFloat(0, fp);

        float pixels[] = (float[]) fp.getPixels();
        for (int i = 0; i < pixels.length; i++) {
            double item = pixels[i] - mean;
            std_dev += item * item;
        }

        std_dev = Math.sqrt(std_dev / pixels.length);

        return std_dev;
    }
}
