
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
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
public class Xmipp_SelReader extends ImagePlus implements PlugIn {

    public void run(String filename) {
        boolean show = false;

        // If launched from menu...
        if (filename == null || filename.isEmpty()) {
            OpenDialog od = new OpenDialog("Load xmipp SEL file...", System.getProperty("user.dir"), "");
            String rootDir = od.getDirectory();
            filename = od.getFileName();

            if (filename == null) {
                return;
            }

            filename = rootDir + filename;
            show = true;
        }

        if (filename == null) {
            return;
        }

        IJ.showStatus("Reading: " + filename);

        try {
            MetaData md = new MetaData();
            md.read(filename);

            if (md.containsLabel(MDLabel.MDL_IMAGE)) {
                long ids[] = md.findObjects();

                ImageStack is = null;

                for (long id : ids) {
                    String imagefilename = md.getValueString(MDLabel.MDL_IMAGE, id);

                    // Reads xmipp image.
                    Xmipp_Reader imp = new Xmipp_Reader();
                    imp.run(imagefilename);

                    if (is == null) {
                        is = new ImageStack(imp.getWidth(), imp.getHeight());
                    }

                    is.addSlice("", imp.getProcessor());
                }

                ImagePlus imp = new ImagePlus(filename, is);

                // Attach the Image Processor
                if (imp.getNSlices() > 1) {
                    setStack(filename, imp.getStack());
                } else {
                    setProcessor(filename, imp.getProcessor());
                }

                // Copy the scale info over
                copyScale(imp);

                // Resets min and max display values.
                getProcessor().resetMinAndMax();

                // Show the image if it was selected by the file
                // chooser, don't if an argument was passed ie
                // some other ImageJ process called the plugin.
                if (show) {
                    show();
                }
            }
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
        }
    }
}
