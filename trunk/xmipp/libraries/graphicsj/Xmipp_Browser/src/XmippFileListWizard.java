
import browser.commandline.Parameters;
import browser.filebrowsers.JDialogXmippFilesList;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippFileListWizard {

    public static void main(String args[]) {
        Parameters parameters = new Parameters(args);

        JDialogXmippFilesList frameBrowser = new JDialogXmippFilesList(
                parameters.directory, parameters);
        frameBrowser.setVisible(true);
    }
}
