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
public class ImageOperations {

    public static ImagePlus mean(Vector<TableImageItem> items) {
        int w = items.get(0).getWidth();
        int h = items.get(0).getHeight();
        double mean[] = new double[w * h];

        // For all images...
        for (int k = 0; k < items.size(); k++) {
            TableImageItem current = items.elementAt(k);

            if (current.isEnabled() && current.exists()) {
                ImagePlus ip = current.getImagePlus();
                float pixels[] = (float[]) ip.getProcessor().getPixels();
                ip.close();

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

    public static ImagePlus std_deviation(Vector<TableImageItem> items) {
        int w = items.get(0).getWidth();
        int h = items.get(0).getHeight();

        float mean[] = (float[]) mean(items).getProcessor().getPixels();
        double stdDev[] = new double[mean.length];

        // Traverses images.
        for (int k = 0; k < items.size(); k++) {
            TableImageItem current = items.elementAt(k);

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
}
