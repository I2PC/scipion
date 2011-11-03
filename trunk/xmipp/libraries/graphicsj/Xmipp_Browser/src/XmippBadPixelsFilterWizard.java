
import browser.commandline.Parameters;
import wizards.JFrameXmippBadPixelsFilter;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippBadPixelsFilterWizard {

    public static void main(String args[]) {
        Parameters parameters = new Parameters(args);

        JFrameXmippBadPixelsFilter frameBrowser = new JFrameXmippBadPixelsFilter(parameters.files[0], parameters);
        frameBrowser.setVisible(true);
    }
}
