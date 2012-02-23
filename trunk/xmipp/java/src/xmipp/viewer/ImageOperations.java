/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer;

import xmipp.viewer.imageitems.tableitems.AbstractGalleryImageItem;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageStatistics;
import java.util.ArrayList;

/**
 *
 * @author Juanjo Vega
 */
public class ImageOperations {

    public static ImagePlus mean(ArrayList<AbstractGalleryImageItem> items) {
        int w = items.get(0).getWidth();
        int h = items.get(0).getHeight();
        double mean[] = new double[w * h];

        // For all images...
        for (AbstractGalleryImageItem current : items) {
            if (current.isEnabled() && current.exists()) {
                ImagePlus ip = current.getImagePlus();
                float pixels[] = (float[]) ip.getProcessor().getPixels();

                // Adds current image to sum.
                for (int i = 0; i < pixels.length; i++) {
                    mean[i] += pixels[i];
                }
            }
        }

        // Calculates mean...
        for (int i = 0; i < mean.length; i++) {
            mean[i] /= items.size(); // ...by dividing with the #images
        }

        return new ImagePlus("Mean", new FloatProcessor(w, h, mean));
    }

    public static ImagePlus std_deviation(ArrayList<AbstractGalleryImageItem> items) {
        int w = items.get(0).getWidth();
        int h = items.get(0).getHeight();

        float mean[] = (float[]) mean(items).getProcessor().getPixels();
        double stdDev[] = new double[mean.length];

        // Traverses images.
        for (AbstractGalleryImageItem current : items) {//int k = 0; k < items.size(); k++) {
            //TableImageItem current = items.elementAt(k);

            if (current.isEnabled() && current.exists()) {
                ImagePlus ip = current.getImagePlus();
                float pixels[] = (float[]) ip.getProcessor().getPixels();
                ip.close();

                for (int i = 0; i < pixels.length; i++) {
                    stdDev[i] += (pixels[i] - mean[i]) * (pixels[i] - mean[i]);
                }
            }
        }

        for (int i = 0; i < stdDev.length; i++) {
            stdDev[i] = Math.sqrt(stdDev[i] / items.size());
        }

        FloatProcessor processor = new FloatProcessor(w, h, stdDev);
        ImagePlus std_ip = new ImagePlus();
        std_ip.setProcessor("Standard Deviation", processor);

        return std_ip;
    }

    public static double[] getMinAndMax(ArrayList<AbstractGalleryImageItem> items) {
        // Initial min and max values.
        ImageStatistics statistics;

        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;

        for (int i = 0; i < items.size(); i++) {
            statistics = items.get(i).getStatistics();

            if (statistics == null) {
                statistics = items.get(i).getPreview().getStatistics();
            }

            if (statistics.min < min) {
                min = statistics.min;
            }

            if (statistics.max > max) {
                max = statistics.max;
            }
        }

        return new double[]{min, max};
    }
    /*public static double[] getMinAndMax(Vector<TableImageItem> items) {
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
    }*/
}
