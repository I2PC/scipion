
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;
import xmipp.Filename;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public abstract class Writer implements PlugIn {

    String basedir, prefix, filename;

    public void run(String arg) {
        ImagePlus imp = WindowManager.getCurrentImage();

        if (imp != null) {
            run(arg, imp);
        } else {
            IJ.error("There are no images open!");
        }
    }

    public void run(String arg, ImagePlus imp) {
        if (imp == null) {
            IJ.error("image is NULL");
            return;
        }

        boolean show = false;
        String path = arg;

        // If was opened by direct call to the plugin
        // instead of through HandleExtraFileTypes, shows dialog.
        if (path == null || path.trim().isEmpty()) {
            path = getPath();
            show = true;
        } else {
            basedir = System.getProperty("user.dir");
            prefix = Filename.getPrefix(path);
            filename = Filename.getFilename(path);
        }

        if (path != null) {
            try {
                write(imp, path);

                if (show) {
                    IJ.showMessage("File sucessfully saved!");
                }
            } catch (Exception ex) {
                ex.printStackTrace();
                IJ.error(ex.getMessage());
            }
        }
    }

//
//    public void run(ImageProcessor ip) {
//        String path = basedir + filename;
//
//        System.out.println("Rootdir: " + basedir);
//        System.out.println("Filename: " + filename);
//
//        IJ.showStatus("Saving: " + filename);
//
//        try {
//            write(imp, path);
//
//            IJ.showMessage("File successfully saved!");
//        } catch (Exception ex) {
//            ex.printStackTrace();
//            IJ.error(ex.getMessage());
//        }
//    }
//
//    public int showDialog(ImagePlus ip, String string, PlugInFilterRunner pifr) {
//        if (filename == null || filename.isEmpty()) {
//            SaveDialog sd = new SaveDialog(getSaveDialogTitle(), System.getProperty("user.dir"), "");
//            rootdir = sd.getDirectory();
//            filename = sd.getFileName();
//            System.out.println("rootdir: " + rootdir + " > " + System.currentTimeMillis());
//            if (filename == null) {
//                return DONE;
//            }
//        }
//
//        return 0;
//    }
    String getPath() {
        String fullpath = null;

        SaveDialog od = new SaveDialog(
                getSaveDialogTitle(), System.getProperty("user.dir"), "", ".xmd");
        String path = od.getFileName();

        if (path != null) {
            basedir = od.getDirectory();
            prefix = Filename.getPrefix(path);
            filename = Filename.getFilename(path);
            System.out.println("basedir : " + basedir);
            System.out.println("prefix  : " + prefix);
            System.out.println("fn      : " + filename);
            fullpath = prefix + Filename.SEPARATOR + basedir + filename;
            System.out.println("fullpath: " + fullpath);
        }

        return fullpath;
    }

    protected abstract String getSaveDialogTitle();

    protected abstract void write(ImagePlus imp, String path) throws Exception;
}
