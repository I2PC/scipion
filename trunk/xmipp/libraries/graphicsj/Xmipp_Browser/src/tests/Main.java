package tests;

import browser.imageitems.ImageConverter;
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
        String path = "/home/juanjo/temp/mask.vol";
        ImageDouble image = new ImageDouble();

        try {
            image.readPreview(path, 80, 80);
            ImageConverter.convertToImagej(image, path).show();
            image.printShape();
            //System.err.println(" *** Nimages: " + image.getNimages());
        } catch (Exception ex) {
            ex.printStackTrace();
            System.err.println(" *** Exception: " + ex.getMessage());
        }
    }
}
