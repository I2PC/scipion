
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
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
public class Xmipp_MetaDataWriter extends Xmipp_Writer {

    @Override
    protected void write(ImagePlus imp, String path) throws Exception {
        // Builds metadata.
        MetaData md = new MetaData();
        md.addLabel(MDLabel.MDL_IMAGE);

        // Calculate name pattern for images.
        String ext = ".xmp";

        String dirname = (new File(path)).getName();
        dirname = dirname.substring(0, dirname.lastIndexOf('.'));
        File dir = new File(rootdir + dirname);
        if (!dir.exists()) {
            dir.mkdir();
        }

        // Saves each slice and adds its name to metadata to save it later.
        ImageStack stack = imp.getStack();
        for (int i = 1; i <= stack.getSize(); i++) {
            String imageName = dir.getAbsolutePath() + File.separator + i + ext;

            ImagePlus slice = new ImagePlus(imageName, stack.getProcessor(i));
            IJ.run(slice, "Xmipp writer", "save=" + imageName);

            // Adds imagename to metadata.
            long id = md.addObject();
            md.setValueString(MDLabel.MDL_IMAGE, imageName, id);
        }

        // Finally, saves the metadata file.
        md.write(path);
    }

    @Override
    protected String getSaveDialogTitle() {
        return "Save as xmipp metadata...";
    }
}
