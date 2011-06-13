
import ij.ImagePlus;
import xmipp.Filename;
import xmipp.ImageDouble;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Xmipp_ImageReader extends Xmipp_Reader {

    @Override
    protected void read(String path) throws Exception {
        ImageDouble image = new ImageDouble(filename);

        ImagePlus imp = ImageConverter.convertToImagej(image, filename);

        // Sets stack...
        String name = Filename.getFilename(filename);
        setStack(name, imp.getStack());

        // ...and copies scale info.
        copyScale(imp);
    }

    @Override
    protected String getOpenDialogTitle() {
        return "Load xmipp file...";
    }
}
