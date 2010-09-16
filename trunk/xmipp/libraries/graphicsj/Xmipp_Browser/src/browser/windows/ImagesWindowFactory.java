/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.windows;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import java.io.File;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesWindowFactory {

    public static ImagePlus openImageWindow(File file) {
        return openImageWindow(file, -1);
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

                openImageWindow(ip2);
            } else {
                // Not an image, so let ImageJ open it as usually.
                IJ.open(file.getAbsolutePath());
            }
        } else {
            IJ.showMessage("Trying to open a null file.");
        }

        return ip;
    }

    public static void openImageWindow(ImagePlus ip) {
        if (ip != null) {
            ImageCanvas cc = new ImageCanvas(ip);

            if (ip.getStackSize() > 1) {
                new WindowStackOperations(ip, cc);
            } else {
                new WindowImageOperations(ip, cc);
            }
        } else {
            IJ.showMessage("Trying to open a null image.");
        }
    }
}
