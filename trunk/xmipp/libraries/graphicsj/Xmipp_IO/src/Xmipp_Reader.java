
import ij.IJ;
import ij.ImagePlus;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import xmipp.ImageDouble;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Xmipp_Reader extends ImagePlus implements PlugIn {

    public void run(String filename) {
        boolean show = false;

        // If launched from menu...
        if (filename == null || filename.isEmpty()) {
            OpenDialog od = new OpenDialog("Load xmipp file...", System.getProperty("user.dir"), "");
            String rootDir = od.getDirectory();
            filename = od.getFileName();

            if (filename == null) {
                return;
            }

            filename = rootDir + filename;
            show = true;
        }

        try {
            ImageDouble image = new ImageDouble(filename);

            ImagePlus imp = ImageConverter.convertToImagej(image, filename);

            // Attach the Image Processor
            if (imp.getNSlices() > 1) {
                setStack(filename, imp.getStack());
            } else {
                setProcessor(filename, imp.getProcessor());
            }

            // Copy the scale info over
            copyScale(imp);

            // Show the image if it was selected by the file
            // chooser, don't if an argument was passed ie
            // some other ImageJ process called the plugin.
            if (show) {
                show();
            }
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
        }
    }
}
