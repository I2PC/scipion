
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileInfo;
import java.io.File;
import xmipp.Filename;
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
public class MetaDataWriter extends Writer {

    @Override
    protected void write(ImagePlus imp, String path) throws Exception {
        FileInfo fi = imp.getOriginalFileInfo();
        File f = new File(fi.directory, fi.fileName);

        String stackfilename = f.getAbsolutePath();
        boolean isMD = Filename.isMetadata(stackfilename);

        if (!f.exists() && !isMD) {
            f = new File(Filename.getFilename(path));
            stackfilename = f.getAbsolutePath();

            stackfilename = stackfilename.substring(0, stackfilename.lastIndexOf("."));
            stackfilename += ".stk";

            ImageWriter iw = new ImageWriter();
            iw.writeStack(imp, stackfilename);
        }

        writeMetaData(imp, stackfilename, isMD, path);
    }

    protected void writeMetaData(ImagePlus imp, String stackFilename, boolean isMD, String path) throws Exception {
        // Builds metadata.
        MetaData md = new MetaData();
        md.addLabel(MDLabel.MDL_IMAGE);

        // Saves each slice and adds its name to metadata to save it later.
        ImageStack stack = imp.getStack();
        for (int i = 1; i <= stack.getSize(); i++) {
            // If MD: getSlice label.
            // Otherwise, builds it.
            String imageName = isMD ? imp.getStack().getSliceLabel(i) : i + Filename.SEPARATOR + stackFilename;

            // Adds imagename to metadata.
            long id = md.addObject();
            md.setValueString(MDLabel.MDL_IMAGE, imageName, id);
        }

        // Finally, saves the metadata file.
        md.write(path);
    }

    @Override
    protected String getSaveDialogTitle() {
        return "Save as MetaData...";
    }
}
