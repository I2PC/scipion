
import ij.IJ;
import ij.ImagePlus;
import ij.io.SaveDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.util.Tools;
import java.io.File;
import xmipp.ImageDouble;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Xmipp_Writer implements PlugInFilter {

    protected String filename, rootdir;
    protected ImagePlus img;

    public int setup(String filename, ImagePlus img) {
        this.filename = filename;
        this.img = img;

        File file = new File(filename);
        rootdir = file.isAbsolute() ? "" : file.getParent() + File.separator;

        return DOES_32 + DOES_16 + DOES_8G + NO_CHANGES;
    }

    public void run(ImageProcessor ip) {
        if (filename == null || filename.isEmpty()) {
            SaveDialog sd = new SaveDialog("save as xmipp...", System.getProperty("user.dir"), "");
            rootdir = sd.getDirectory();
            filename = sd.getFileName();
            if (filename == null) {
                return;
            }
        }

        IJ.showStatus("Saving: " + rootdir + filename);

        try {
            float data[] = (float[]) img.getProcessor().getPixels();
            write(filename,
                    img.getWidth(), img.getHeight(), img.getStackSize(),
                    Tools.toDouble(data));
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
        }
    }

    protected static void write(String filename, int w, int h, int d, double data[]) throws Exception {
        ImageDouble image = new ImageDouble();

        image.setFilename(filename);
        image.setData(w, h, d, data);
        image.write(filename);
    }
}
