/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.windows;

import ij.IJ;
import ij.ImagePlus;
import java.io.File;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesWindowFactory {

    public static ImagePlus openImageWindow(File file) {
        return openImageWindow(file, -1);
    }

    public static ImagePlus openImageWindow(String file, boolean poll) {
        File f = new File(file);
        ImagePlus ip = IJ.openImage(f.getAbsolutePath());
        openImageWindow(ip, poll);

        return ip;
    }

    public static ImagePlus openImageWindow(File file, int slice) {
        ImagePlus ip = null;

        if (file != null) {
            ip = IJ.openImage(file.getAbsolutePath());

            if (ip != null) {
                ImagePlus ip2 = ip;

                if (slice >= 0) {
                    ip2 = new ImagePlus();
                    ip2.getProcessor().setPixels(ip.getStack().getPixels(slice));
                }

                openImageWindow(ip2, false);
            } else {
                // Not an image, so let ImageJ open it as usually.
                IJ.open(file.getAbsolutePath());
            }
        } else {
            IJ.showMessage("Trying to open a null file.");
        }

        return ip;
    }

    public static void openImageWindow(ImagePlus ip, boolean poll) {
        if (ip != null) {
            if (ip.getStackSize() > 1) {
                new WindowStackOperations(ip, poll);
            } else {
                new WindowImageOperations(ip, poll);
            }
        } else {
            IJ.showMessage("Trying to open a null image.");
        }
    }
}
