
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

    public static boolean saveResult(List<Patch> patches, Rectangle box, int margin, String path) {

        try {
            // Translation to the top left corner.
            AffineTransform translation = new AffineTransform();
            translation.translate(box.x, box.y);

            // Creates a new ImagePlus big enough.
            int w = box.width;
            int h = box.height;

            float result[] = new float[w * h];

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

            Rectangle area0 = new Rectangle(image0.getWidth(), image0.getHeight());
            Rectangle area1 = new Rectangle(image1.getWidth(), image1.getHeight());

            // Transforms concatenated with translation to the box top left corner.
            AffineTransform T0 = patch0.getAffineTransform().createInverse();
            AffineTransform T1 = patch1.getAffineTransform().createInverse();

            T0.concatenate(translation);
            T1.concatenate(translation);

            Point p = new Point();
            Point p0 = new Point();
            Point p1 = new Point();

            // Gets statistics to normalize later.
            System.out.println(" >>> Retrieving overlapped area mean and stdDev.");
            float meanAndStdDev[] = getStatistics(
                    pixels0, width0, height0, area0, T0,
                    pixels1, width1, height1, area1, T1,
                    w, h, margin);
            float mean0 = meanAndStdDev[0];
            float stdDev0 = meanAndStdDev[1];
            float mean1 = meanAndStdDev[2];
            float stdDev1 = meanAndStdDev[3];

            // Assigns value.
            System.out.println(" >>> Building result image.");
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    p.setLocation(x, y);

                    // Transforms both points: p -> p1, p2
                    T0.transform(p, p0);
                    T1.transform(p, p1);

                    float value = -1;
                    float value0 = 0;   // TODO: Replace
                    float value1 = 1;   // TODO: Replace
                    if (area0.contains(p0) && area1.contains(p1)) {   // Overlap!
                        // Horizontal blending.
                        float w0 = (float) (width0 - p0.getX());
                        float w1 = (float) p1.getX();

                        // Vertical blending.
                        float h0 = (float) p0.getY();
                        float h1 = (float) (height1 - p1.getY());

                        w0 = Math.abs(w0);
                        w1 = Math.abs(w1);
                        h0 = Math.abs(h0);
                        h1 = Math.abs(h1);

                        /*float value0 = Interpolator.bilinear(pixels0, width0, height0, p0.getX(), p0.getY());
                        float value1 = Interpolator.bilinear(pixels1, width1, height1, p1.getX(), p1.getY());
                        
                        // Normalizes values.
                        value0 = (float) ((value0 - mean0) / stdDev0);
                        value1 = (float) ((value1 - mean1) / stdDev1);*/

                        // Blending.
                        if (w0 > margin && h0 > margin) {
                            value = value0;
                        } else {
                            // Horizontal blending.
                            float w0_ = w0 / (w0 + w1);
                            float w1_ = w1 / (w0 + w1);

                            // Vertical blending.
                            float h0_ = h0 / (h0 + h1);
                            float h1_ = h1 / (h0 + h1);

                            // Combined blending.
                            w0_ = (w0_ * h0_) / (w0_ + h1_);
                            w1_ = (w1_ * h1_) / (w0_ + h1_);

                            // Smoothstep blending: t*t*(3 - 2*t)
                            w0_ = w0_ * w0_ * (3 - 2 * w0_);
                            w1_ = w1_ * w1_ * (3 - 2 * w1_);

                            value = w0_ * value0 + w1_ * value1;
                        }
                    } else if (area0.contains(p0)) {   // Contained only by Image0.
                        //value = Interpolator.bilinear(pixels0, width0, height0, p0.getX(), p0.getY());
                        //value = (float) ((value - mean0) / stdDev0);
                        value = value0;
                    } else if (area1.contains(p1)) {   // Contained only by Image1.
                        //value = Interpolator.bilinear(pixels1, width1, height1, p1.getX(), p1.getY());
                        //value = (float) ((value - mean1) / stdDev1);
                        value = value1;
                    }

                    // Normalizes value.
                    result[y * w + x] = value;
                }
            }

            // Writes slice to disk.            
            System.err.println(" >>> Saving result: " + path);
            ImageDouble image = new ImageDouble();
            image.setData(box.width, box.height, 1, toDouble(result));
            image.write(path);

            new ImageJ();
            ImageDouble image2 = new ImageDouble(path);
            double data[] = image2.getData();
            FloatProcessor fp = new FloatProcessor(box.width, box.height, data);
            ImagePlus imp = new ImagePlus(path, fp);
            imp.show();
        } catch (Exception ex) {
            ex.printStackTrace(System.err);
            IJ.error("Exception: " + ex.getMessage());

            return false;
        }

        return true;
    }

    static float[] getStatistics(float pixels0[], int width0, int height0, Rectangle area0, AffineTransform T0,
            float pixels1[], int width1, int height1, Rectangle area1, AffineTransform T1, int w, int h, float margin) {
        Point p = new Point();
        Point p0 = new Point();
        Point p1 = new Point();

        float mean0 = 0, mean1 = 0;
        float stdDev0 = 0, stdDev1 = 0;
        int n = 0;

        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                p.setLocation(x, y);

                // Transforms both points: p -> p1, p2
                T0.transform(p, p0);
                T1.transform(p, p1);

                if (area0.contains(p0) && area1.contains(p1)) {   // Overlap!
                    float value0 = Interpolator.bilinear(pixels0, width0, height0, p0.getX(), p0.getY());
                    float value1 = Interpolator.bilinear(pixels1, width1, height1, p1.getX(), p1.getY());

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
        stdDev0 = (float) Math.sqrt(stdDev0);
        stdDev1 = stdDev1 / n - mean1 * mean1;
        stdDev1 = (float) Math.sqrt(stdDev1);

        return new float[]{mean0, stdDev0, mean1, stdDev1};
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

        Point p = new Point();

        for (int j = 0; j < h; j++) {   // Rows
            p.y = j;
            for (int i = 0; i < w; i++) {   // Columns
                p.x = i;

                Point p_ = new Point();
                T.transform(p, p_);

                // If point is inside image...
                if (r.contains(p_)) {

                    // Bilinear Interpolation.
                    float value = Interpolator.bilinear(pixels, ip.getWidth(), ip.getHeight(), p_.getX(), p_.getY());

                    // Store point.
                    slice[j * w + i] = value;
                }
            }
        }

        return slice;
    }
}
