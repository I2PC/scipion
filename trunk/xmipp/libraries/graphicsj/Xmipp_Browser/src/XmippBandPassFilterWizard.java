
import browser.commandline.Parameters;
import wizards.JFrameXmippBandPassFilter;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippBandPassFilterWizard {

    public static void main(String args[]) {
        Parameters parameters = new Parameters(args);

        JFrameXmippBandPassFilter frameBrowser = new JFrameXmippBandPassFilter(
                parameters.files[0], parameters);
        frameBrowser.setVisible(true);
    }
}
