/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tests;

import browser.windows.ImagesWindowFactory;
import ij.IJ;
import ij.ImageJ;
import ij.plugin.PlugIn;
import java.io.File;
import javax.swing.JFileChooser;

/**
 *
 * @author Juanjo Vega
 */
public class PollTest implements PlugIn {

    public PollTest() {
    }

    public static void main(String args[]) {
        new ImageJ();

//        (new PollTest()).run("");
        ImagesWindowFactory.openImageWindow(IJ.openImage("/home/juanjo/Desktop/Stack.tif"), false);
    }

    public void run(String string) {
        JFileChooser jfc = new JFileChooser();

        if (jfc.showOpenDialog(null) != JFileChooser.CANCEL_OPTION) {
            File f = jfc.getSelectedFile();

            // Opens image using poll (reloads continuously from disk)
            ImagesWindowFactory.openImageWindow(IJ.openImage(f.getAbsolutePath()), true);
        }
    }
}
