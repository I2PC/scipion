
import ij.ImagePlus;
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
public class ImageWriter extends Writer {

    @Override
    protected void write(ImagePlus imp, String path) throws Exception {
        // Stack is fully written.       
        String fullPath = imp.getStackSize() > 1 ? Filename.getFilename(path) : path;
        writeStack(imp, fullPath);
    }

    protected void writeStack(ImagePlus imp, String path) throws Exception {
        int w = imp.getWidth();
        int h = imp.getHeight();
        int d = imp.getStackSize();

        for (int i = 0; i < d; i++) {
            float slice[] = (float[]) imp.getStack().getProcessor(i + 1).getPixels();

            ImageDouble image = new ImageDouble();
            image.setData(w, h, 1, toDouble(slice));
            image.write(i + Filename.SEPARATOR + path);
        }
    }

    double[] toDouble(float input[]) {
        double output[] = new double[input.length];

        for (int i = 0; i < input.length; i++) {
            output[i] = input[i];
        }

        return output;
    }

    @Override
    protected String getSaveDialogTitle() {
        return "Save as xmipp image...";
    }
}
