
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileInfo;
import java.io.File;
import xmipp.MetaData;

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
    protected void read(String path) throws Exception {
        IJ.showStatus("Reading: " + path);

        MetaData md = new MetaData(path);

        ImagePlus imp = ImageConverter.convertToImagej(md);

        File f = new File(path);
        FileInfo fi = new FileInfo();
        fi.directory = f.getParent();
        fi.fileName = f.getName();
        setFileInfo(fi);

        // Sets stack...
        setStack(getTitle(), imp.getStack());

        // ...and copies scale info.
        copyScale(imp);
    }

    @Override
    protected String getOpenDialogTitle() {
        return "Open MetaData...";
    }
}
