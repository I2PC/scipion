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

    public static ImagePlus mean(Vector<TableImageItem> items) {
        ImagePlus current = items.firstElement().getImagePlus();
        int w = current.getWidth();
        int h = current.getHeight();
        float sum[][] = new float[w][h];

        // For all images...
        for (int k = 0; k < items.size(); k++) {
            current = items.elementAt(k).getImagePlus();

            // Adds current image to sum.
            for (int j = 0; j < h; j++) {
                for (int i = 0; i < w; i++) {
                    sum[i][j] += current.getProcessor().getPixelValue(i, j);
                }
            }

            current.close();
            current = items.elementAt(k).getImagePlus();
        }

        float mean[][] = new float[w][h];

        // Calculates mean...
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                mean[i][j] = sum[i][j] / items.size(); // ...by dividing with the #images
            }
        }

        FloatProcessor processor = new FloatProcessor(mean);
        ImagePlus meanIp = new ImagePlus();
        meanIp.setProcessor("", processor);

        return meanIp;
    }

    public static ImagePlus std_deviation(Vector<TableImageItem> items) {
        ImagePlus current = items.firstElement().getImagePlus();
        int w = current.getWidth();
        int h = current.getHeight();

        ImagePlus mean = mean(items);
        float std[] = new float[mean.getWidth() * mean.getHeight()];

        // Traverses images.
        for (int k = 0; k < items.size(); k++) {
            current = items.elementAt(k).getImagePlus();

            for (int j = 0; j < h; j++) {
                for (int i = 0; i < w; i++) {
                    float pixel = current.getProcessor().getPixelValue(i, j);
                    float avg = mean.getProcessor().getPixelValue(i, j);

                    std[j * w + i] += (pixel - avg) * (pixel - avg);
                }
            }
            current.close();
        }

        for (int i = 0; i < h; i++) {
            for (int j = 0; j < h; j++) {
                std[j * w + i] /= items.size();
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
