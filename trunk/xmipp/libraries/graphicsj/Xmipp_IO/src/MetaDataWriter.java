
import ij.IJ;
import ij.ImageJ;
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

    public static void main(String args[]) {
        new ImageJ();

        String filename = "class_000002@/data2/MicrographPreprocessing/2D/CL2D/run_001/results_level_00_classes.xmd";

        try {
            MetaDataReader xr = new MetaDataReader();
            xr.run(filename);
            xr.show();

            String output = "/home/juanjo/Desktop/stackfull.xmd";
            MetaDataWriter iw = new MetaDataWriter();
            iw.write(xr, output);

            // Reload.
//            ImageReader ir = new ImageReader();
//            ir.run(output);
//            ir.show();
        } catch (Exception ex) {
            ex.printStackTrace();
            IJ.error(ex.getMessage());
        }
    }

    @Override
    protected void write(ImagePlus imp, String path) throws Exception {
        FileInfo fi = imp.getFileInfo();
        String directory = fi.directory;
        String filename = fi.fileName;

        System.out.println(" > directory: " + directory);
        System.out.println(" > fileName: " + filename);

        // @TODO @TODO @TODO @TODO @TODO @TODO
        // If file !exists
        // Write as stack: path(-ext)+".stk"
        String stackfilename = "";
        File f = new File(directory, filename);
        if (!f.exists()) {
            stackfilename = f.getPath();

            stackfilename = stackfilename.substring(0, stackfilename.lastIndexOf("."));
            stackfilename += ".stk";

            ImageWriter iw = new ImageWriter();
            iw.writeStack(imp, stackfilename);
        }

        writeMetaData(imp, stackfilename, path);
    }

    protected void writeMetaData(ImagePlus imp, String stackFilename, String path) throws Exception {
        // Builds metadata.
        MetaData md = new MetaData();
        md.addLabel(MDLabel.MDL_IMAGE);

        // Saves each slice and adds its name to metadata to save it later.
        ImageStack stack = imp.getStack();
        for (int i = 1; i <= stack.getSize(); i++) {
            String imageName = i + Filename.SEPARATOR + stackFilename;

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
