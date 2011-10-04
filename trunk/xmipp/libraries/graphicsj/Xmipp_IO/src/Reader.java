
import ij.IJ;
import ij.ImagePlus;
import ij.io.OpenDialog;
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
public abstract class Reader extends ImagePlus implements PlugIn {

    String basedir, prefix, filename;

    public void run(String arg) {
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
        String fullpath = null;

        OpenDialog od = new OpenDialog(getOpenDialogTitle(), System.getProperty("user.dir"), "");
        String path = od.getFileName();

        if (path != null) {
            basedir = od.getDirectory();
            prefix = Filename.getPrefix(path);
            filename = Filename.getFilename(path);
            fullpath = prefix + Filename.SEPARATOR + basedir + filename;
        }

        return fullpath;
    }

    @Override
    public String getTitle() {
        String p = prefix == null ? "" : prefix + Filename.SEPARATOR;

        return p + filename;
    }

    protected abstract String getOpenDialogTitle();

    protected abstract void read(String path) throws Exception;
}
