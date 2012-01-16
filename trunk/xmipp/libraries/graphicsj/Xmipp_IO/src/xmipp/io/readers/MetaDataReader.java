package xmipp.io.readers;

import ij.IJ;
import ij.ImagePlus;
import ij.io.FileInfo;
import java.io.File;
import xmipp.jni.MetaData;
import xmippij.XmippImageConverter;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class MetaDataReader extends Reader {

    @Override
    public void read(String path) throws Exception {
        IJ.showStatus("Reading: " + path);

        MetaData md = new MetaData(path);

        ImagePlus imp = XmippImageConverter.convertToImageJ(md);

        File f = new File(path);
        FileInfo fi = new FileInfo();
        fi.directory = f.getParent();
        fi.fileName = f.getName();
        setFileInfo(fi);

        // Sets stack...
        setStack(imp.getTitle(), imp.getStack());

        // ...and copies scale info.
        copyScale(imp);
    }

    @Override
    protected String getOpenDialogTitle() {
        return "Open MetaData...";
    }
}
