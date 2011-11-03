
import browser.commandline.Parameters;
import wizards.JFrameXmippFilesListMeasureFreqs;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippFileListMeasureFreqsWizard {

    public static void main(String args[]) {
        Parameters parameters = new Parameters(args);

        JFrameXmippFilesListMeasureFreqs frameBrowser = new JFrameXmippFilesListMeasureFreqs(
                parameters.directory, parameters);
        frameBrowser.setVisible(true);
    }
}
