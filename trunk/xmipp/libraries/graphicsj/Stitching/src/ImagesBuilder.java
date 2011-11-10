
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ini.trakem2.display.Patch;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.List;
import xmipp.Filename;
import xmipp.ImageDouble;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class ImagesBuilder {

    static boolean DEBUG = true;

    public static boolean saveResult(List<Patch> patches, Rectangle box, int margin, String path) {

        try {
            // Translation to the top left corner.
            AffineTransform translation = AffineTransform.getTranslateInstance(box.x, box.y);

            // Creates a new ImagePlus big enough.
            int w = box.width;
            int h = box.height;

            double result[] = new double[w * h];

            // Puts images, pixels, translations in arrays, a much convenient way to use them later.
            Patch patch0 = patches.get(0);
            Patch patch1 = patches.get(1);

            ImagePlus image0 = patch0.getImagePlus();
            ImagePlus image1 = patch1.getImagePlus();

            int width0 = image0.getWidth();
            int height0 = image0.getHeight();
            int width1 = image1.getWidth();
            int height1 = image1.getHeight();

            float pixels0[] = (float[]) image0.getProcessor().convertToFloat().getPixels();
            float pixels1[] = (float[]) image1.getProcessor().convertToFloat().getPixels();

            Rectangle area0 = new Rectangle(width0, height0);
            Rectangle area1 = new Rectangle(width1, height1);

            // Upper-left and bottom-right corners.
            Point p00 = new Point(0, 0);
            Point p01 = new Point(width0, height0);
            Point p10 = new Point(0, 0);
            Point p11 = new Point(width1, height1);

            Point.Double p = new Point.Double();
            Point.Double p0 = new Point.Double();
            Point.Double p1 = new Point.Double();

            // Translations for images points.
            AffineTransform T0 = patch0.getAffineTransformCopy().createInverse();
            AffineTransform T1 = patch1.getAffineTransformCopy().createInverse();
            T0.concatenate(translation);
            T1.concatenate(translation);

            // Gets statistics to normalize later.
            System.out.println(" >>> Retrieving overlapped area mean and stdDev.");
            double meanAndStdDev[] = getStatistics(
                    pixels0, width0, height0, area0, T0,
                    pixels1, width1, height1, area1, T1,
                    w, h/*, margin*/);
            double mean0 = meanAndStdDev[0];
            double stdDev0 = meanAndStdDev[1];
            double mean1 = meanAndStdDev[2];
            double stdDev1 = meanAndStdDev[3];

            // Blending.
            double wx, wy;
            double value0, value1;

            // Assigns value.
            System.out.println(" >>> Building result image.");
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    p.setLocation(x, y);

                    // Transforms both points: p -> p1, p2
                    T0.transform(p, p0);
                    T1.transform(p, p1);

                    double d0x = Math.min(p0.getX(), Math.abs(p01.getX() - p0.getX()));
                    double d0y = Math.min(p0.getY(), Math.abs(p01.getY() - p0.getY()));

                    value0 = value1 = 0;
                    wx = wy = 0;

                    if (area0.contains(p0)) {
                        value0 = Interpolator.bilinear(pixels0, width0, height0, p0.getX(), p0.getY());

                        if (area1.contains(p1)) {   // Overlap!
                            value1 = Interpolator.bilinear(pixels1, width1, height1, p1.getX(), p1.getY());

                            if (d0x > margin) {
                                wx = 1;
                            } else {
                                wx = d0x / margin;
                            }

                            if (d0y > margin) {
                                wy = 1;
                            } else {
                                wy = d0y / margin;
                            }
                        } else {    // Image 0
                            wx = wy = 1;
                        }
                    } else if (area1.contains(p1)) {    // Image 1
                        value1 = Interpolator.bilinear(pixels1, width1, height1, p1.getX(), p1.getY());
                        wx = wy = 0;
                    }

                    // Normalizes values.
                    value0 = (value0 - mean0) / stdDev0;
                    value1 = (value1 - mean1) / stdDev1;

                    double w0 = wx * wy;
                    double w1 = 1 - wx * wy;
                    double isum = 1 / (w0 + w1);

                    // Smoothstep.
//                    wx = wx * wx * (3 - 2 * wx);
//                    wy = wy * wy * (3 - 2 * wy);
                    wx = wx * wx * wx * (wx * (wx * 6 - 15) + 10);
                    wy = wy * wy * wy * (wy * (wy * 6 - 15) + 10);

                    result[y * w + x] = isum * (w0 * value0 + w1 * value1);
                }
            }

            // Writes slice to disk.            
            System.out.println(" >>> Saving result: " + path);
            FloatProcessor fp = new FloatProcessor(w, h, result);
            ImagePlus imp = new ImagePlus("Result", fp);
            IJ.save(imp, path);

//            new ImageJ();
//            IJ.open(path);
        } catch (Exception ex) {
            ex.printStackTrace(System.err);
            IJ.error("Exception: " + ex.getMessage());

            return false;
        }

        return true;
    }

    /**
     * Returns the minimum distance between a point and two lines
     * |   |
     * a p b
     * |   |
     * @param p point
     * @param a line a
     * @param b line b (opposite to a)
     * @return minimum between distance(p, a) and distance(p, b)
     */
    /*    static double getMinimumDistance(Point.Double p, Line2D a, Line2D b) {
    double da = a.ptLineDist(p);
    double db = b.ptLineDist(p);
    return Math.min(da, db);
    }*/
    static double[] getStatistics(float pixels0[], int width0, int height0, Rectangle area0, AffineTransform T0,
            float pixels1[], int width1, int height1, Rectangle area1, AffineTransform T1, int w, int h/*, double margin*/) {
        Point p = new Point();
        Point p0 = new Point();
        Point p1 = new Point();

        double mean0 = 0, mean1 = 0;
        double stdDev0 = 0, stdDev1 = 0;
        int n = 0;

        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                p.setLocation(x, y);

                // Transforms both points: p -> p1, p2
                T0.transform(p, p0);
                T1.transform(p, p1);

                if (area0.contains(p0) && area1.contains(p1)) {   // Overlap!
                    double value0 = Interpolator.bilinear(pixels0, width0, height0, p0.getX(), p0.getY());
                    double value1 = Interpolator.bilinear(pixels1, width1, height1, p1.getX(), p1.getY());

                    mean0 += value0;
                    mean1 += value1;

                    stdDev0 += value0 * value0;
                    stdDev1 += value1 * value1;

                    n++;
                }
            }
        }

        // Mean and Standard deviation.
        mean0 /= n;
        mean1 /= n;

        stdDev0 = stdDev0 / n - mean0 * mean0;
        stdDev0 = (double) Math.sqrt(stdDev0);
        stdDev1 = stdDev1 / n - mean1 * mean1;
        stdDev1 = (double) Math.sqrt(stdDev1);

        return new double[]{mean0, stdDev0, mean1, stdDev1};
    }

    public static boolean saveStack(List<Patch> patches, Rectangle box, String path) {
        File f = new File(path);

        // Removes previous (to avoid adding to a stack with different size images)
        if (f.exists()) {
            f.delete();
        }

        for (int i = 1; i <= patches.size(); i++) {
            Patch patch = patches.get(i - 1);

            try {
                float slice[] = buildSlice(patch, box);

                ImageDouble image = new ImageDouble();
                image.setData(box.width, box.height, 1, toDouble(slice));

                // Writes slice to disk.
                System.err.println(" >>> Saving slice: " + i + Filename.SEPARATOR + path);
                image.write(i + Filename.SEPARATOR + path);
            } catch (Exception ex) {
                ex.printStackTrace(System.err);
                IJ.error("Exception: " + ex.getMessage());
                return false;
            }
        }

        return true;
    }

    static double[] toDouble(float[] input) {
        double output[] = null;

        if (input != null) {
            output = new double[input.length];
            for (int i = 0; i < input.length; i++) {
                output[i] = input[i];
            }
        }

        return output;
    }

    static float[] buildSlice(Patch patch, Rectangle box) throws Exception {
        ImagePlus ip = patch.getImagePlus();
        AffineTransform translation = new AffineTransform();
        translation.translate(box.x, box.y);

        AffineTransform impT = patch.getAffineTransform();
        impT.invert(); // Inverts transformation.

        AffineTransform T = new AffineTransform(impT);
        T.concatenate(translation);

        float pixels[] = (float[]) ip.getProcessor().convertToFloat().getPixels();
        Rectangle r = new Rectangle(ip.getWidth(), ip.getHeight());

        BufferedImage bimage = new BufferedImage(box.width, box.height,
                BufferedImage.TYPE_INT_ARGB);
        ImagePlus imp = new ImagePlus("tmp", bimage);

        int w = box.width;
        int h = box.height;
        float slice[] = (float[]) imp.getProcessor().convertToFloat().getPixels();

        Point.Double p = new Point.Double();

        for (int j = 0; j < h; j++) {   // Rows
            p.y = j;
            for (int i = 0; i < w; i++) {   // Columns
                p.x = i;

                Point.Double p_ = new Point.Double();
                T.transform(p, p_);

                // If point is inside image...
                if (r.contains(p_)) {

                    // Bilinear Interpolation.
                    double value = Interpolator.bilinear(pixels, ip.getWidth(), ip.getHeight(), p_.getX(), p_.getY());

                    // Store point.
                    slice[j * w + i] = (float) value;
                }
            }
        }

        return slice;
    }
}
