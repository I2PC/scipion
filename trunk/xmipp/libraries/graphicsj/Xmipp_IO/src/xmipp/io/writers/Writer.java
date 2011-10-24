package xmipp.io.writers;

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
                ex.printStackTrace(System.err);
                IJ.error(ex.getMessage());
            }
        }
    }

    String getPath() {
        String fullpath = null;

        SaveDialog od = new SaveDialog(
                getSaveDialogTitle(), System.getProperty("user.dir"), "", "");
        String path = od.getFileName();

        if (path != null) {
            basedir = od.getDirectory();
            prefix = Filename.getPrefix(path);
            prefix = prefix != null ? prefix + Filename.SEPARATOR : "";
            filename = Filename.getFilename(path);
            fullpath = prefix + basedir + filename;
        }

        return fullpath;
    }

    protected abstract String getSaveDialogTitle();

    public abstract void write(ImagePlus imp, String path) throws Exception;
}
