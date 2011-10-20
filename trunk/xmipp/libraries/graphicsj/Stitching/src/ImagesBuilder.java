
import ij.IJ;
import ij.ImagePlus;
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
    /*
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
    }*/

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

            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    p.setLocation(x, y);

                    // Transforms both points: p -> p1, p2
                    T0.transform(p, p0);
                    T1.transform(p, p1);

                    /*                    System.out.println(p + " > " + p0);
                    System.out.println(p + " > " + p1);
                    System.out.println("- - - - - - - - - - - - -");*/

                    float value = 0;
                    if (area0.contains(p0) && area1.contains(p1)) {   // Overlap!
                        float w0 = (float) (width0 - p0.getX());
                        float w1 = (float) p1.getX();

                        float value0 = Interpolator.bilinear(pixels0, width0, height0, p0.getX(), p0.getY());
                        float value1 = Interpolator.bilinear(pixels1, width1, height1, p1.getX(), p1.getY());

                        // TODO Usar una propiedad en el archivo para el "margin".
                        if (w0 > margin) {
                            value = value0;
                        } else {
                            float w0_ = w0 / (w0 + w1);
                            float w1_ = w1 / (w0 + w1);

                            // "Blends" value.
                            value = w0_ * value0 + w1_ * value1;
                        }
                    } else if (area0.contains(p0)) {   // Contained only by Image0
                        value = Interpolator.bilinear(pixels0, width0, height0, p0.getX(), p0.getY());
                    } else if (area1.contains(p1)) {   // Contained only  by Image1
                        value = Interpolator.bilinear(pixels1, width1, height1, p1.getX(), p1.getY());
                    }

                    result[y * w + x] = value;
                }
            }

            // Writes slice to disk.            
            System.err.println(" >>> Saving result: " + path);
            ImageDouble image = new ImageDouble();
            image.setData(box.width, box.height, 1, toDouble(result));
            image.write(path);
        } catch (Exception ex) {
            ex.printStackTrace(System.err);
            IJ.error("Exception: " + ex.getMessage());

            return false;
        }

        return true;
    }

    /**
     * Blend a FloatProcessor (i.e., a 32-bit image) with another one, i.e.
     * set the pixel values of fp1 according to a weighted sum of the corresponding
     * pixels of fp1 and fp2. This is done for pixels in the rectangle fp1.getRoi()
     * only.
     * Note that both FloatProcessors, fp1 and fp2 must have the same width and height.
     * @param fp1 The FloatProcessor that will be modified.
     * @param weight1 The weight of the pixels of fp1 in the sum.
     * @param fp2  The FloatProcessor that will be read only.
     * @param weight2 The weight of the pixels of fp2 in the sum.
     */
//    public static void blendFloat (FloatProcessor fp1, float weight1, FloatProcessor fp2, float weight2) {
//        Rectangle r = fp1.getRoi();
//        int width = fp1.getWidth();
//        float[] pixels1 = (float[])fp1.getPixels();     // array of the pixels of fp1
//        float[] pixels2 = (float[])fp2.getPixels();
//        for (int y=r.y; y<(r.y+r.height); y++)          // loop over all pixels inside the roi rectangle
//            for (int x=r.x; x<(r.x+r.width); x++) {
//                int i = x + y*width;                    // this is how the pixels are addressed
//                pixels1[i] = weight1*pixels1[i] + weight2*pixels2[i]; //the weighted sum
//            }
//    }

    /*static int getRectangleIndexContaining(Rectangle rectangles[], Point point) {
    int contains = -1;
    int i = 0;
    
    for (; i < rectangles.length; i++) {
    Rectangle rectangle = rectangles[i];
    if (rectangle.contains(point)) {
    if (contains >= 0) {    // More than one rectangle contains the point.
    return -1;
    }
    
    contains = i;
    }
    }
    
    return contains;
    }*/
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
//System.err.println(" -> X: " + image.getXsize()+ " --> "+box.width);
//System.err.println(" -> Y: " + image.getYsize()+ " --> "+box.height);
//System.err.println(" -> Z: " + image.getNsize()+ " / "+patches.size());

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
    /*
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
    }*/
}
