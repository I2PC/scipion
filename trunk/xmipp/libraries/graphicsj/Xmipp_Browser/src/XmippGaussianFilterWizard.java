
import browser.commandline.Parameters;
import wizards.JFrameXmippGaussianFilter;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippGaussianFilterWizard {

    public static void main(String args[]) {
        Parameters parameters = new Parameters(args);

        JFrameXmippGaussianFilter frameBrowser = new JFrameXmippGaussianFilter(
                parameters.files[0], parameters);
        frameBrowser.setVisible(true);
    }
}
