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
        double average[] = null;
        int w = 0, h = 0;

        // For all images...
        for (int k = 0; k < items.size(); k++) {
            TableImageItem current = items.elementAt(k);

            if (current.isEnabled()) {
                ImagePlus ip = current.getImagePlus();

                if (average == null) {
                    w = current.getWidth();
                    h = current.getHeight();
                    average = new double[w * h];
                }

                // Adds current image to sum.
                for (int j = 0; j < h; j++) {
                    for (int i = 0; i < w; i++) {
                        average[j * w + i] += ip.getProcessor().getPixelValue(i, j);
                    }
                }

                ip.close();
            }
        }

        // Calculates average...
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                average[j * w + i] /= items.size(); // ...by dividing with the #images
            }
        }

        FloatProcessor processor = new FloatProcessor(w, h, average);
        ImagePlus averageIp = new ImagePlus();
        averageIp.setProcessor("Average", processor);

        return averageIp;
    }

    public static ImagePlus std_deviation(Vector<TableImageItem> items) {
        double stdDev[] = null;
        int w = 0, h = 0;

        ImagePlus average = average(items);

        // Traverses images.
        for (int k = 0; k < items.size(); k++) {
            TableImageItem current = items.elementAt(k);

            if (current.isEnabled()) {
                ImagePlus ip = current.getImagePlus();

                if (stdDev == null) {
                    w = current.getWidth();
                    h = current.getHeight();
                    stdDev = new double[w * h];
                }

                for (int j = 0; j < h; j++) {
                    for (int i = 0; i < w; i++) {
                        float pixel = ip.getProcessor().getPixelValue(i, j);
                        float avg = average.getProcessor().getPixelValue(i, j);

                        stdDev[j * w + i] += (pixel - avg) * (pixel - avg);
                    }
                }
                ip.close();
            }
        }

        for (int i = 0; i < h; i++) {
            for (int j = 0; j < h; j++) {
                stdDev[j * w + i] = Math.sqrt(stdDev[j * w + i] / items.size());
            }
        }

        FloatProcessor processor = new FloatProcessor(w, h, stdDev);
        ImagePlus std_ip = new ImagePlus();
        std_ip.setProcessor("Std. Dev.", processor);

        return std_ip;
    }

    public static double[] getMinAndMax(Vector<TableImageItem> items) {
        // Initial min and max values.
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

    public static double averagePixelValue(ImagePlus imp) {
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

    public static double stdDevPixelValue(ImagePlus imp) {
        double mean = averagePixelValue(imp);
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
