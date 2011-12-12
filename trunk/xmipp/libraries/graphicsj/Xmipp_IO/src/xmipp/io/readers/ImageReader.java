package xmipp.io.readers;

import ij.IJ;
import ij.ImagePlus;
import ij.io.FileInfo;
import java.io.File;
import xmipp.Filename;
import xmippij.XmippImageConverter;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class ImageReader extends Reader {

    @Override
    public void read(String path) throws Exception {
        IJ.showStatus("Reading: " + path);
        System.out.println("Reading: " + path);

/*        String fullpath = path;
        if (prefix != null && !prefix.trim().isEmpty()) {
            long n = Long.valueOf(prefix);
            fullpath = n + Filename.SEPARATOR + path;
        }

        System.out.println(" *** fullpath: " + fullpath);
*/
        ImagePlus imp = XmippImageConverter.loadImage(path);

        // @TODO Extract path, dir, etc.
        File f = new File(path);
        FileInfo fi = new FileInfo();
        fi.directory = f.getParent();
        fi.fileName = f.getName();
        setFileInfo(fi);

        // Sets stack...
        String name = Filename.getFilename(path);
        setStack(name, imp.getStack());

        // ...and copies scale info.
        copyScale(imp);
    }

    @Override
    protected String getOpenDialogTitle() {
        return "Open xmipp image...";
    }
}
