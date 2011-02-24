package sphere;


/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Interpolator {

    public static double bilinear(SphereCell table[][], double x, double y) {
        int x_ = (int) x;
        int y_ = (int) y;

        double v00 = table[x_][y_].value;
        double v01 = table[x_][y_].value;
        double v10 = table[x_][y_].value;
        double v11 = table[x_][y_].value;

        return bilinear(v00, v01, v10, v11, x, y);
    }

    /**
     * Bilinear interpolation
     * @param v00 i,j value
     * @param v01 i+1,j value
     * @param v10 i+1,j value
     * @param v11 i+1,j+1 value
     * @param x point to interpolate
     * @return
     */
    private static double bilinear(double v00, double v01, double v10, double v11, double x, double y) {
        double v_up = linear(v00, v01, x);
        double v_down = linear(v10, v11, x);

        return linear(v_up, v_down, y);
    }

    /**
     * Linear interpolation
     * @param v_left i value
     * @param v_right i+1 value
     * @param x point to interpolate
     * @return
     */
    protected static double linear(double v_left, double v_right, double x) {
        int int_x = (int) x;

        double d = x - int_x;

        return v_right * d - v_left * (d - 1);
    }
}
