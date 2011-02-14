/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tests;

import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;
import javax.swing.JToggleButton;

/**
 *
 * @author Juanjo Vega
 */
public class ImageNormalizationTest implements PlugIn {

    public static void main(String args[]) {
        new ImageJ();
        (new ImageNormalizationTest()).run("");
    }

    public void run(String string) {
        String dir = "/home/juanjo/Desktop/selfiles/images";
        String fileNames[] = {
            "12205_000001.xmp",
            "12205_000002.xmp",
            "12205_000003.xmp",
            "12205_000004.xmp",
            "12205_000005.xmp"};

        final Vector<ImagePlus> images = new Vector<ImagePlus>();
        final JToggleButton b = new JToggleButton("Normalize");

        b.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                if (b.isSelected()) {

                    double parameters[] = getMultipleNormalizationParameters(images);

                    double min = parameters[0];
                    double max = parameters[1];
                    double factor = (double) 255 / (double) (max - min);

                    for (int j = 0; j < images.size(); j++) {
                        System.out.println(
                                j + ": " + images.elementAt(j).getTitle()
                                + " / m=" + images.elementAt(j).getProcessor().getMin()
                                + " / M=" + images.elementAt(j).getProcessor().getMax());

                        images.elementAt(j).getProcessor().setMinAndMax(min, max);
                        images.elementAt(j).updateAndDraw();

                        System.out.println(
                                j + ": " + images.elementAt(j).getTitle()
                                + " / m=" + images.elementAt(j).getProcessor().getMin()
                                + " / M=" + images.elementAt(j).getProcessor().getMax());
                    }
                } else {
                    for (int j = 0; j < images.size(); j++) {
                        images.elementAt(j).getProcessor().resetMinAndMax();
                        images.elementAt(j).updateAndDraw();
                    }
                }
            }
        });

        for (int i = 0; i < fileNames.length; i++) {
//            Spider_Reader sr = new Spider_Reader();
            ImagePlus image = /*sr.load(dir, fileNames[i]);  */ createIP(i, i + 1);

            image.show();

            if (i == fileNames.length - 1) {
                image.getWindow().add(b);
                image.getWindow().pack();
            }

            images.add(image);
        }
    }

    private ImagePlus createIP(float min, float max) {
        int w = 100, h = 100;
        ImageProcessor ip = new FloatProcessor(w, h);
        final ImagePlus image = new ImagePlus("Test", ip);

        float[] pixels = (float[]) ip.getPixels();

        // Create custom image.
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                pixels[j * w + i] = (j < h / 2) ? min : max;
            }
        }

        ip.setPixels(pixels);

        return image;
    }

    public static double[] getMultipleNormalizationParameters(Vector<ImagePlus> items) {
        // Calculates min and max.
        double min = items.elementAt(0).getProcessor().getMin();
        double max = items.elementAt(0).getProcessor().getMax();

        for (int i = 1; i < items.size(); i++) {
            double current_min = items.elementAt(i).getProcessor().getMin();
            double current_max = items.elementAt(i).getProcessor().getMax();

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

        return new double[]{min, max, factor};
    }
}
