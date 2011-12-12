package xmipp.io.readers;

import ij.IJ;
import ij.ImagePlus;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public abstract class Reader extends ImagePlus implements PlugIn {

    String path = "";

    public void run(String arg) {
        boolean show = false;
        path = arg;

        // If was opened by direct call to the plugin
        // instead of through HandleExtraFileTypes, shows dialog.
        if (path == null || path.trim().isEmpty()) {
            path = getPath();
            show = true;
        }

        if (path != null) {
            try {
                read(path);

                if (show) {
                    show();
                }
            } catch (Exception ex) {
                ex.printStackTrace(System.err);
                IJ.error(ex.getMessage());
            }
        }
    }

    String getPath() {
        OpenDialog od = new OpenDialog(getOpenDialogTitle(), System.getProperty("user.dir"), "");

        return od.getFileName();
    }

    @Override
    public String getTitle() {
        return path;
    }

    protected abstract String getOpenDialogTitle();

    public abstract void read(String path) throws Exception;
}
