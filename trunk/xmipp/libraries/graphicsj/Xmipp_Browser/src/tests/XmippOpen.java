package tests;

import browser.imageitems.ImageConverter;
import ij.IJ;
import java.io.File;
import xmipp.ImageDouble;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippOpen {

    public static void main(String args[]) {
        //new ImageJ();

        (new XmippOpen()).run(args);
    }

    public void run(String args[]) {
        try {
            String path = "/media/PENDRIVE/Ad5GLflagIIIa/a09251/down3_a09251.spi";//args[0];  // "/home/juanjo/Desktop/imgs_Roberto/g7106.raw";
            if (!path.startsWith(File.separator)) {
                path = System.getProperty("user.dir") + File.separator + path;
            }

            System.out.println(" *** Loading: [" + path + "]");

            ImageDouble image = new ImageDouble();
            image.readHeader(path);

            int w = image.getWidth() / 4;
            int h = image.getHeight() / 4;
            System.out.println(" *** w=" + w + " h=" + h);

            image = new ImageDouble();
            image.readPreview(path, w, h, 1);
            //image.read(path);

            ImageConverter.convertToImagej(image, path).show();
            System.out.println(" > Done!");
        } catch (Exception ex) {
            ex.printStackTrace();
            IJ.error(ex.getMessage());
        }
    }
}
