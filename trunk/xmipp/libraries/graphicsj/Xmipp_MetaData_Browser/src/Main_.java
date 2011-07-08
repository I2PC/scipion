
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import browser.DEBUG;
import ij.ImageJ;
import metadata.JFrameMetaData;

/**
 *
 * @author Juanjo Vega
 */
public class Main_ {

//    static {
//        System.setProperty("plugins.dir", "/home/juanjo/Desktop/ImageJ/plugins");
//    }
    public static void main(String args[]) {
        try {
            new ImageJ();

            runMetaDataBrowser(args);
        } catch (Exception ex) {
            ex.printStackTrace();
            //throw new RuntimeException(ex);
        }
    }

    public static void runMetaDataBrowser(String args[]) {
        String filename = "/data2/coss/Preprocessing/micrographs_blocks.sel";

        DEBUG.printMessage("Loading micrograph: " + filename);

        JFrameMetaData frameMetaData = new JFrameMetaData(filename);
        frameMetaData.setLocationRelativeTo(null);
        frameMetaData.setVisible(true);
    }
}
