
import ij.IJ;
import ij.ImagePlus;
import ij.io.SaveDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import java.io.File;
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
public class Xmipp_SelWriter implements PlugInFilter {

    protected String filename, rootdir;
    protected ImagePlus img;

    public int setup(String filename, ImagePlus img) {
        this.filename = filename;
        this.img = img;

        File file = new File(filename);
        rootdir = file.getParent() + File.separator;

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

        System.out.println("Saving: " + rootdir + filename);

        IJ.showMessage("@TODO Save!! + ABOUT");

        try {
            // Builds metadata.
            MetaData md = new MetaData();
            md.addLabel(MDLabel.MDL_IMAGE);

            // @TODO Convert reversed into ImageDouble: for EACH slice.
            for (int i = 0; i < img.getStackSize(); i++) {
                String imageName = filename + i;
                // ImagePlus imp = ImageConverter.convertToImagej(image, filename);
                ImageDouble image = new ImageDouble();

                image.write(imageName); // Writes image to disk.

                // Add imagename to metadata.
                long id = md.addObject();
                md.setValueString(MDLabel.MDL_IMAGE, imageName, id);
            }

            // Finally, writes the metadata file.
            md.write(filename);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
        }
    }
}
