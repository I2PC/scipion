
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.SaveDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.util.Tools;
import java.io.File;
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
public class Xmipp_SelWriter implements PlugInFilter {

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
            SaveDialog sd = new SaveDialog("save as xmipp SEL...", System.getProperty("user.dir"), "");
            rootdir = sd.getDirectory();
            filename = sd.getFileName();

            if (filename == null) {
                return;
            }
        }

        IJ.showStatus("Saving: " + rootdir + filename);

        try {
            // Builds metadata.
            MetaData md = new MetaData();
            md.addLabel(MDLabel.MDL_IMAGE);

            // Calculate name pattern for images.
            String name = filename.substring(0, filename.lastIndexOf('.'));
            String ext = ".xmp";

            // Saves each slice and adds its name to metadata to save it later.
            ImageStack stack = img.getStack();
            for (int i = 1; i <= img.getStackSize(); i++) {
                String imageName = name + i + ext;

                float data[] = (float[]) stack.getProcessor(i).getPixels();
                Xmipp_Writer.write(imageName,
                        img.getWidth(), img.getHeight(), 1,
                        Tools.toDouble(data));
                //run("Xmipp writer", "save=/home/juanjo/Desktop/kk.xmp");

                // Add imagename to metadata.
                long id = md.addObject();
                md.setValueString(MDLabel.MDL_IMAGE, imageName, id);
            }

            // Finally, saves the metadata file.
            md.write(filename);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
        }
    }
}
