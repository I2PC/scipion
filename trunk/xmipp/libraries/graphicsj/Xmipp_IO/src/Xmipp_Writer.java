
import ij.IJ;
import ij.ImagePlus;
import ij.io.SaveDialog;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
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
public abstract class Xmipp_Writer implements ExtendedPlugInFilter {

    protected String filename, rootdir;
    protected ImagePlus imp;

    public int setup(String filename, ImagePlus imp) {
        this.filename = filename;
        this.imp = imp;

        File file = new File(filename);
        rootdir = file.isAbsolute() ? "" : file.getParent() + File.separator;

        return DOES_32 + DOES_16 + DOES_8G + NO_CHANGES;
    }

    public int showDialog(ImagePlus ip, String string, PlugInFilterRunner pifr) {
        if (filename == null || filename.isEmpty()) {
            SaveDialog sd = new SaveDialog(getSaveDialogTitle(), System.getProperty("user.dir"), "");
            rootdir = sd.getDirectory();
            filename = sd.getFileName();

            if (filename == null) {
                return DONE;
            }
        }

        return 0;
    }

    public void setNPasses(int i) {
    }

    public void run(ImageProcessor ip) {
        String path = rootdir + filename;

        IJ.showStatus("Saving: " + filename);

        try {
            write(imp, path);

            IJ.showMessage("File successfully saved!");
        } catch (Exception ex) {
            ex.printStackTrace();
            IJ.error(ex.getMessage());
        }
    }

    protected abstract String getSaveDialogTitle();

    protected abstract void write(ImagePlus imp, String path) throws Exception;
}
