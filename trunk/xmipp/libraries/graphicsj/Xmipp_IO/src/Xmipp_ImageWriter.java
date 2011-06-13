
import ij.ImagePlus;
import ij.util.Tools;
import xmipp.ImageDouble;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Xmipp_ImageWriter extends Xmipp_Writer {

    @Override
    protected void write(ImagePlus imp, String path) throws Exception {
        int w = imp.getWidth();
        int h = imp.getHeight();
        int d = imp.getStackSize();

        ImageDouble image = new ImageDouble();
        image.setFilename(path);

        // Copies the entire stack.
        double data[] = new double[w * h * d];
        for (int i = 0; i < d; i++) {
            float slice[] = (float[]) imp.getStack().getProcessor(i + 1).getPixels();
            System.arraycopy(Tools.toDouble(slice), 0, data, i * w * h, slice.length);
        }

        image.setData(w, h, d, data);
        image.write(path);
    }

    @Override
    protected String getSaveDialogTitle() {
        return "save as xmipp...";
    }
}
