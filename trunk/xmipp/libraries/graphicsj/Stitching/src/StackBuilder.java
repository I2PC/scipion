
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ini.trakem2.display.Patch;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.util.List;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class StackBuilder {

    public static ImagePlus buildStack(List<Patch> patches, Rectangle box) {
        ImageStack is = new ImageStack(box.width, box.height);

        for (int i = 0; i < patches.size(); i++) {
            Patch patch = patches.get(i);

            try {
                float slice[] = buildSlice(patch, box);

                is.addSlice(patch.getImagePlus().getTitle(), slice);
            } catch (Exception e) {
                IJ.log("Exception: " + e.getMessage());
            }
        }

        ImagePlus ip = new ImagePlus("", is);

        return ip;
    }

    private static float[] buildSlice(Patch patch, Rectangle box) throws Exception {
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
        float slice[] = (float[]) imp.getProcessor().convertToFloat().getPixels();//new float[w * h];

        Point p = new Point();

        for (int j = 0; j < h; j++) {
            p.y = j;
            for (int i = 0; i < w; i++) {
                p.x = i;

                Point p_ = new Point();
                T.transform(p, p_);

                // If point is inside image...
                if (r.contains(p_)) {

                    // Bilinear Interpolation.
                    float value = Interpolator.bilinear(pixels, ip.getWidth(), ip.getHeight(), p_.getX(), p_.getY());

                    // Store point.
                    slice[j * w + i] = value;//pixels[y * ip.getWidth() + x];
                }
            }
        }

        return slice;
    }

    public static ImagePlus sumStack(ImagePlus imp) {
        int w = imp.getWidth();
        int h = imp.getHeight();
        float pixels[] = new float[w * h];

        for (int i = 0; i < imp.getStackSize(); i++) {
            float current[] = (float[]) imp.getStack().getProcessor(i + 1).getPixels();

            pixels = sumArrays(pixels, current);
        }

        FloatProcessor ip = new FloatProcessor(w, h);
        ip.setPixels(pixels);

        return new ImagePlus("", ip);
    }

    static float[] sumArrays(float[] a, float[] b) {
        float result[] = new float[a.length];

        for (int i = 0; i < result.length; i++) {
            result[i] = a[i] + b[i];
        }

        return result;
    }
}
