
import browser.commandline.Parameters;
import wizards.JFrameXmippMaskDesign;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippMaskDesignWizard {

    public static void main(String args[]) {
        Parameters parameters = new Parameters(args);

        JFrameXmippMaskDesign frameBrowser = new JFrameXmippMaskDesign(
                parameters.files[0], parameters);
        frameBrowser.setVisible(true);
    }
}
