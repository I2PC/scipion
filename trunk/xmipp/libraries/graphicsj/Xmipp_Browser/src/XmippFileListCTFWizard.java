
import browser.commandline.Parameters;
import wizards.JFrameXmippFilesListCTF;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippFileListCTFWizard {

    public static void main(String args[]) {
        Parameters parameters = new Parameters(args);

        JFrameXmippFilesListCTF frameBrowser = new JFrameXmippFilesListCTF(
                parameters.directory, parameters);
        frameBrowser.setVisible(true);
    }
}
