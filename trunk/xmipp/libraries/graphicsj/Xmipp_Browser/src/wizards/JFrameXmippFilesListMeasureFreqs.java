/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.LABELS;
import browser.commandline.Parameters;
import browser.filebrowsers.JPanelXmippBrowser;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameXmippFilesListMeasureFreqs extends JFrameXmippFilesListCTF {

    public JFrameXmippFilesListMeasureFreqs(String directory, Parameters parameters) {
        super(directory, parameters);

        setTitle(LABELS.TITLE_WIZARD_MEASURE_FREQS);

        jpButtons.remove(jbOk);
        jbCancel.setText(LABELS.BUTTON_CLOSE);
    }

    @Override
    protected JPanelXmippBrowser createPanel(Parameters parameters) {
        JPanelXmippFileListCTF panel = new JPanelXmippFileListCTF(parameters.directory,
                parameters.filter, parameters.downsampling);
        panel.setDownsamplingEnabled(false);

        return panel;
    }
}
