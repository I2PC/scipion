
import browser.filebrowsers.JFrameXmippFilesListPSD;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippFileListPSD extends XmippFileList {

    void runBrowser(String directory, String expression, boolean singleSelection) {
        JFrameXmippFilesListPSD frameBrowser = new JFrameXmippFilesListPSD(
                directory, PORT, expression, singleSelection);
        frameBrowser.setVisible(true);
    }
}
