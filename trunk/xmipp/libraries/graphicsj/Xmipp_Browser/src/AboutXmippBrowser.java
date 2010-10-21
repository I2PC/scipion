
import browser.LABELS;
import ij.IJ;
import ij.plugin.PlugIn;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class AboutXmippBrowser implements PlugIn {

    public void run(String string) {

        try {
            IJ.showStatus("Starting \"" + LABELS.APP_NAME + "\" plugin...");
            IJ.showMessage("About " + LABELS.APP_NAME, "@TO DO");
        } catch (RuntimeException e) {
            String msg = LABELS.APP_NAME + "\n" +
                    "*****************************************************************\n" +
                    "Original error message:\n" + e;
            IJ.showMessage(LABELS.APP_NAME, msg);
        } finally {
            IJ.showStatus("");
        }
    }
}
