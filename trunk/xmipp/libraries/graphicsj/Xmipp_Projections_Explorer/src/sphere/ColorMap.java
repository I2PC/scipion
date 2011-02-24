package sphere;


import sphere.Sphere;
import ij.ImagePlus;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import java.awt.Color;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class ColorMap {

    // Color map Jet: 64 colors palette.
    private final static double MAP[][] = {
        {0, 0, 0.562500000000000},
        {0, 0, 0.625000000000000},
        {0, 0, 0.687500000000000},
        {0, 0, 0.750000000000000},
        {0, 0, 0.812500000000000},
        {0, 0, 0.875000000000000},
        {0, 0, 0.937500000000000},
        {0, 0, 1},
        {0, 0.0625000000000000, 1},
        {0, 0.125000000000000, 1},
        {0, 0.187500000000000, 1},
        {0, 0.250000000000000, 1},
        {0, 0.312500000000000, 1},
        {0, 0.375000000000000, 1},
        {0, 0.437500000000000, 1},
        {0, 0.500000000000000, 1},
        {0, 0.562500000000000, 1},
        {0, 0.625000000000000, 1},
        {0, 0.687500000000000, 1},
        {0, 0.750000000000000, 1},
        {0, 0.812500000000000, 1},
        {0, 0.875000000000000, 1},
        {0, 0.937500000000000, 1},
        {0, 1, 1},
        {0.0625000000000000, 1, 0.937500000000000},
        {0.125000000000000, 1, 0.875000000000000},
        {0.187500000000000, 1, 0.812500000000000},
        {0.250000000000000, 1, 0.750000000000000},
        {0.312500000000000, 1, 0.687500000000000},
        {0.375000000000000, 1, 0.625000000000000},
        {0.437500000000000, 1, 0.562500000000000},
        {0.500000000000000, 1, 0.500000000000000},
        {0.562500000000000, 1, 0.437500000000000},
        {0.625000000000000, 1, 0.375000000000000},
        {0.687500000000000, 1, 0.312500000000000},
        {0.750000000000000, 1, 0.250000000000000},
        {0.812500000000000, 1, 0.187500000000000},
        {0.875000000000000, 1, 0.125000000000000},
        {0.937500000000000, 1, 0.0625000000000000},
        {1, 1, 0},
        {1, 0.937500000000000, 0},
        {1, 0.875000000000000, 0},
        {1, 0.812500000000000, 0},
        {1, 0.750000000000000, 0},
        {1, 0.687500000000000, 0},
        {1, 0.625000000000000, 0},
        {1, 0.562500000000000, 0},
        {1, 0.500000000000000, 0},
        {1, 0.437500000000000, 0},
        {1, 0.375000000000000, 0},
        {1, 0.312500000000000, 0},
        {1, 0.250000000000000, 0},
        {1, 0.187500000000000, 0},
        {1, 0.125000000000000, 0},
        {1, 0.0625000000000000, 0},
        {1, 0, 0},
        {0.937500000000000, 0, 0},
        {0.875000000000000, 0, 0},
        {0.812500000000000, 0, 0},
        {0.750000000000000, 0, 0},
        {0.687500000000000, 0, 0},
        {0.625000000000000, 0, 0},
        {0.562500000000000, 0, 0},
        {0.500000000000000, 0, 0}};

    public static int getColor(double index, double max) {
        double i = index * (MAP.length / max);  // scales index because map has only 64 values.

        int i_int = (int) i;

        // Left and right colors.
        double left[] = MAP[i_int];
        double right[] = i_int < MAP.length - 1 ? MAP[i_int + 1] : MAP[i_int];

        // Interpolator for each R G B pair.
        double v[] = new double[3];
        for (int j = 0; j < v.length; j++) {
            v[j] = Sphere.VALUE * Interpolator.linear(left[j], right[j], i);
        }

//        System.out.println(i + " > " + i_int + ": " + v[0] / Sphere.VALUE + ", " + v[1] / Sphere.VALUE + ", " + v[2] / Sphere.VALUE);

        // Finally builds color.
        return (new Color(
                (int) Math.round(v[0]),
                (int) Math.round(v[1]),
                (int) Math.round(v[2]))).getRGB();
    }

    /**
     * Builds an ImagePlus to show the map
     */
    public static void show() {
        int H = 256;    // 256 colors.
        int W = 50;     // width
        int map[] = ColorMap.getMap(W, H);

        ImageProcessor ip = new ColorProcessor(W, H);
        ImagePlus image = new ImagePlus("Color Map", ip);
        image.getProcessor().setPixels(map);
        image.show();
    }

    private static int[] getMap(int W, int H) {
        int img[] = new int[H * W];

        for (int i = 0; i < H; i++) {
            int color = getColor(i, H);

            // Copies color into the entire row.
            for (int j = 0; j < W; j++) {
                img[i * W + j] = color;
            }
        }

        return img;
    }
}
