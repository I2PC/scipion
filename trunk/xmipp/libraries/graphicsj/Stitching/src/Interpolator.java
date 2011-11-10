
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Interpolator {

    public static double bilinear(float pixels[], int w, int h, double x, double y) {
        int x_ = (int) x;
        int y_ = (int) y;

        // @TODO Fix bounds
        double v00 = (x_ >= 0 && x_ < w && y_ >= 0 && y_ < h) ? pixels[y_ * w + x_] : 0;
        double v01 = (x_ + 1 >= 0 && x_ + 1 < w && y_ >= 0 && y_ < h) ? pixels[y_ * w + x_ + 1] : 0;
        double v10 = (x_ >= 0 && x_ < w && y_ + 1 >= 0 && y_ + 1 < h) ? pixels[(y_ + 1) * w + x_] : 0;
        double v11 = (x_ + 1 >= 0 && x_ + 1 < w && y_ + 1 >= 0 && y_ + 1 < h) ? pixels[(y_ + 1) * w + (x_ + 1)] : 0;

        return bilinear(v00, v01, v10, v11, (double) x, (double) y);
    }

    /**
     * Bilinear interpolation
     * @param v00 i,j value
     * @param v01 i+1,j value
     * @param v10 i+1,j value
     * @param v11 i+1,j+1 value
     * @param x,y point to interpolate
     * @return
     */
    private static double bilinear(double v00, double v01, double v10, double v11, double x, double y) {
        double v0 = linear(v00, v01, x);
        double v1 = linear(v10, v11, x);

        return linear(v0, v1, y);
    }

    /**
     * Linear interpolation
     * @param v0 i value
     * @param v1 i+1 value
     * @param x value to interpolate
     * @return
     */
    protected static double linear(double v0, double v1, double x) {
        int int_x = (int) x;

        double d = x - int_x;

        return v1 * d - v0 * (d - 1);
    }
}
