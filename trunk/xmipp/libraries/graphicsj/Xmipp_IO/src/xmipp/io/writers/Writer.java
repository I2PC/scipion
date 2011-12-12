package xmipp.io.writers;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public abstract class Writer implements PlugIn {

    public void run(String arg) {
        ImagePlus imp = WindowManager.getCurrentImage();

        if (imp != null) {
            run(imp, arg);
        } else {
            IJ.error("There are no images open!");
        }
    }

    public void run(ImagePlus imp, String arg) {
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
        }

        if (path != null) {
            try {
                write(imp, path);

                imp.changes = false;    // Avoid the "Save changes?" dialog.
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
        SaveDialog od = new SaveDialog(
                getSaveDialogTitle(), System.getProperty("user.dir"), "", "");
        return od.getFileName();
    }

    protected abstract String getSaveDialogTitle();

    public abstract void write(ImagePlus imp, String path) throws Exception;
}
