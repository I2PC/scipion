
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Interpolator {

    public static float bilinear(float pixels[], int w, int h, double x, double y) {
        int x_ = (int) x;
        int y_ = (int) y;

        // @TODO Fix bounds
        float v00 = (x_ >= 0 && x_ < w && y_ >= 0 && y_ < h) ? pixels[y_ * w + x_] : 0;
        float v01 = (x_ + 1 >= 0 && x_ + 1 < w && y_ >= 0 && y_ < h) ? pixels[y_ * w + x_ + 1] : 0;
        float v10 = (x_ >= 0 && x_ < w && y_ + 1 >= 0 && y_ + 1 < h) ? pixels[(y_ + 1) * w + x_] : 0;
        float v11 = (x_ + 1 >= 0 && x_ + 1 < w && y_ + 1 >= 0 && y_ + 1 < h) ? pixels[(y_ + 1) * w + (x_ + 1)] : 0;

        return bilinear(v00, v01, v10, v11, (float) x, (float) y);
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
    private static float bilinear(float v00, float v01, float v10, float v11, float x, float y) {
        float v0 = linear(v00, v01, x);
        float v1 = linear(v10, v11, x);

        return linear(v0, v1, y);
    }

    /**
     * Linear interpolation
     * @param v0 i value
     * @param v1 i+1 value
     * @param x value to interpolate
     * @return
     */
    protected static float linear(float v0, float v1, float x) {
        int int_x = (int) x;

        float d = x - int_x;

        return v1 * d - v0 * (d - 1);
    }
}
