/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.LABELS;
import browser.commandline.Parameters;
import browser.filebrowsers.JDialogXmippFilesList;
import browser.filebrowsers.JPanelXmippBrowser;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameXmippMaskDesign extends JDialogXmippFilesList {

    public JFrameXmippMaskDesign(String metadata, Parameters parameters) {
        super(metadata, parameters);

        setTitle(LABELS.TITLE_WIZARD_MASK_DESIGN);
    }

    @Override
    protected JPanelXmippBrowser createPanel(Parameters parameters) {
        return new JPanelXmippMaskDesign(parameters.files[0], parameters.maskFilename);
    }

    @Override
    protected boolean sendSelectedFiles() {
        return true;
    }
}
