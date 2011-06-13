
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.StackStatistics;
import java.util.LinkedList;
import xmipp.ImageDouble;
import xmipp.MDLabel;
import xmipp.MetaData;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class ImageConverter {

    public static ImagePlus convertToImagej(ImageDouble image, String path) {
        return convertToImagej(image, path, true);
    }

    public static ImagePlus convertToImagej(ImageDouble image, String path, boolean useLogarithm) {
        if (image.isPSD()) {
            image.convertPSD(useLogarithm);
        }

        int w = image.getXsize();
        int h = image.getYsize();
        int d = image.getZsize();
        long n = image.getNsize();

        return convertToImagej(image.getData(), w, h, d, n, path);
    }

    private static ImagePlus convertToImagej(double array[], int w, int h, int d, long n, String title) {
        int sliceSize = w * h;
        int imageSize = sliceSize * d;
        double out[] = new double[sliceSize];
        ImageStack is = new ImageStack(w, h);

        for (int i = 0; i < n; i++) {
            int offset = i * imageSize;
            for (int j = 0; j < d; j++) {
                System.arraycopy(array, offset + j * sliceSize, out, 0, sliceSize);

                FloatProcessor processor = new FloatProcessor(w, h, out);
                is.addSlice(String.valueOf(i), processor);
            }
        }

        ImagePlus ip = new ImagePlus(title, is);

        // Normalize by default
        StackStatistics ss = new StackStatistics(ip);
        ip.getProcessor().setMinAndMax(ss.min, ss.max);

        return ip;
    }

    public static ImagePlus convertToImagej(MetaData md) {
        LinkedList<String> missing = new LinkedList<String>();
        ImagePlus ip = null;

        if (md.containsLabel(MDLabel.MDL_IMAGE)) {
            ImageStack is = null;

            long ids[] = md.findObjects();

            for (long id : ids) {
                String filename = md.getValueString(MDLabel.MDL_IMAGE, id);

                try {
                    ImageDouble img = new ImageDouble(filename);
                    ImagePlus slice = ImageConverter.convertToImagej(img, filename);

                    if (is == null) {
                        is = new ImageStack(slice.getWidth(), slice.getHeight());
                    }

                    is.addSlice(filename, slice.getProcessor());
                } catch (Exception ex) {
                    missing.add(filename);
                    //ex.printStackTrace();
                }
            }

            ip = new ImagePlus(md.getFilename(), is);

            // Normalize by default
            StackStatistics ss = new StackStatistics(ip);
            ip.getProcessor().setMinAndMax(ss.min, ss.max);
        }

        // Tells user about missing files.
        if (!missing.isEmpty()) {
            String message = "There are missing files:\n";
            for (int i = 0; i < missing.size(); i++) {
                message += missing.get(i) + "\n";
            }

            IJ.error(message);
        }

        return ip;
    }
}
