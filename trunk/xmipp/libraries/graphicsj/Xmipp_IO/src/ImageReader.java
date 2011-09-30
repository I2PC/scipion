
import ij.IJ;
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
public class ImageReader extends Reader {

    @Override
    protected void read(String path) throws Exception {
        IJ.showStatus("Reading: " + path);

        ImageDouble image = new ImageDouble();
        if (prefix != null && !prefix.trim().isEmpty()) {
            long n = Long.valueOf(prefix);
            image.read(path, n);
        } else {
            image.read(path);
        }

        ImagePlus imp = ImageConverter.convertToImagej(image, getTitle());

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
