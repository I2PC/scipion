
import ij.IJ;
import ij.ImagePlus;
import xmipp.Filename;
import xmipp.MetaData;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Xmipp_MetaDataReader extends Xmipp_Reader {

    @Override
    protected void read(String path) throws Exception {
        IJ.showStatus("Reading: " + filename);

        try {
            MetaData md = new MetaData();
            md.read(filename);

            ImagePlus imp = ImageConverter.convertToImagej(md);

            // Sets stack...
            String name = Filename.getFilename(filename);
            setStack(name, imp.getStack());

            // ...and copies scale info.
            copyScale(imp);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            ex.printStackTrace();
        }
    }

    @Override
    protected String getOpenDialogTitle() {
        return "Load xmipp MetaData file...";
    }
}
