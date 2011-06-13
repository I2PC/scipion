
import ij.IJ;
import ij.ImagePlus;
import ij.io.OpenDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import java.io.File;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public abstract class Xmipp_Reader extends ImagePlus implements PlugInFilter {

    protected String filename, rootdir;

    public int setup(String filename, ImagePlus ip) {
        this.filename = filename;

        File file = new File(filename);
        rootdir = file.isAbsolute() ? "" : file.getParent() + File.separator;

        return NO_IMAGE_REQUIRED + NO_CHANGES;
    }

    public int showDialog(ImagePlus ip, String string) {
        if (filename == null || filename.isEmpty()) {
            OpenDialog od = new OpenDialog(getOpenDialogTitle(), System.getProperty("user.dir"), "");
            rootdir = od.getDirectory();
            filename = od.getFileName();

            if (filename == null) {
                return -1;
            }
        }

        return 0;
    }

    public void run(ImageProcessor ip) {
        String path = rootdir + filename;

        IJ.showStatus("Reading: " + filename);

        try {
            read(path);
        } catch (Exception ex) {
            ex.printStackTrace();
            IJ.error(ex.getMessage());
        }
    }

    protected abstract String getOpenDialogTitle();

    protected abstract void read(String path) throws Exception;
}
