package xmipp.io.writers;

import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileInfo;
import java.io.File;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class MetaDataWriter extends Writer {
/*
    public static void main(String args[]) {
        try {
            ImagePlus imp = buildMD();
            String output = "/home/jvega/Escritorio/out.xmd";

            MetaDataWriter mdw = new MetaDataWriter();
            mdw.write(imp, output);

            MetaDataReader mdr = new MetaDataReader();
            mdr.run(output);
            mdr.show();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public static ImagePlus buildMD() throws Exception {
        String dir = "/home/jvega/Escritorio/imgs_Roberto/kk/";
        String filenames[] = new String[]{
            dir + "img000001.xmp",
            dir + "img000002.xmp",
            dir + "img000003.xmp",
            dir + "img000004.xmp",
            dir + "img000005.xmp",};

        MetaData md = new MetaData();
        md.addLabel(MDLabel.MDL_IMAGE);

        for (int i = 0; i < filenames.length; i++) {
            String filename = filenames[i];
            long id = md.addObject();

            md.setValueString(MDLabel.MDL_IMAGE, filename, id);
        }

        return XmippImageConverter.convertToImageJ(md);
    }*/

    @Override
    public void write(ImagePlus imp, String path) throws Exception {
        FileInfo fi = imp.getOriginalFileInfo();
        File f = new File(fi.directory, fi.fileName);

        String stackfilename = f.getAbsolutePath();
        boolean isMD = Filename.isMetadata(stackfilename);

        boolean writeStack = fi.fileName == null || fi.fileName.isEmpty() || (!f.exists() && !isMD);

        if (writeStack) {
            f = new File(Filename.getFilename(path));
            stackfilename = f.getAbsolutePath();

            stackfilename = stackfilename.substring(0, stackfilename.lastIndexOf("."));
            stackfilename += ".stk";

            // Writes the entire stack.
            ImageWriter iw = new ImageWriter();
            for (long i = ImageGeneric.FIRST_IMAGE; i <= imp.getStackSize(); i++) {
                ImagePlus slice = new ImagePlus("", imp.getStack().getProcessor((int) i));
                iw.write(slice, i + Filename.SEPARATOR + stackfilename);
            }
        }

        writeMetaData(imp, stackfilename, isMD, path);
    }

    protected void writeMetaData(ImagePlus imp, String stackFilename, boolean isMD, String path) throws Exception {
        // Builds metadata.
        MetaData md = new MetaData();
        md.addLabel(MDLabel.MDL_IMAGE);

        // Saves each slice and adds its name to metadata to save it later.
        ImageStack stack = imp.getStack();
        for (int i = ImageGeneric.FIRST_SLICE; i <= stack.getSize(); i++) {
            // If MD: getSlice label.
            // Otherwise, builds it.
            String imageName = isMD ? stack.getSliceLabel(i) : i + Filename.SEPARATOR + stackFilename;

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
