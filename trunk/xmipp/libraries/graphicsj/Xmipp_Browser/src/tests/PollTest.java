/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tests;

import browser.windows.ImagesWindowFactory;
import ij.IJ;
import ij.ImageJ;
import ij.plugin.PlugIn;

/**
 *
 * @author Juanjo Vega
 */
public class PollTest implements PlugIn {

    public static void main(String args[]) {
        new ImageJ();

        (new PollTest()).run("");
    }

    public void run(String string) {
        // Opens image using poll (reloads continuously from disk)
        ImagesWindowFactory.openImageWindow(IJ.openImage("/home/juanjo/Desktop/poll_test.png"), true);
    }
}
