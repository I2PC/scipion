package tests;

import browser.Cache;
import browser.imageitems.ImageConverter;
import browser.imageitems.listitems.XmippImageItem;
import ij.ImageJ;
import ij.ImagePlus;
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
public class Main {

    public static void main(String args[]) {
        new ImageJ();
        Cache cache = new Cache();

        String path = "/media/PENDRIVE/SampleData/1GRLfullCorrected_050.vol";//"/home/juanjo/temp/Iter_6_filtered_reconstruction.vol";
        //"/media/PENDRIVE/Ad5GLflagIIIa/a09250/down3_a09250_Periodogramavg_enhanced.xmp";
        ImageDouble imageheader = new ImageDouble();

        try {
            imageheader.readHeader(path);
            imageheader.printShape();

            XmippImageItem xii = new XmippImageItem(new File(path), cache);

            xii.getPreview(80, 80).show();
/*            for (int i = 1; i <= 2; i++) {
                ImageDouble image = new ImageDouble();
                image.readPreview(path, 101, 101, i);
                //image.read(path);
//                image.printShape();

                ImagePlus ip = ImageConverter.convertToImagej(image, path);
                ip.setTitle("" + i);
                ip.show();
            }*/
        } catch (Exception ex) {
            ex.printStackTrace();
            System.err.println(" *** Exception: " + ex.getMessage());
        }
    }
}
