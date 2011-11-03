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
public class JFrameXmippBandPassFilter extends JDialogXmippFilesList {

    public JFrameXmippBandPassFilter(String metadata, Parameters parameters) {
        super(metadata, parameters);

        setTitle(LABELS.TITLE_WIZARD_BAND_PASS_FILTER);
//
//        // Hack: Replaces panel.
//        remove(panelXmippBrowser);
//
//        panelXmippBrowser = new JPanelXmippBandPassFilter(metadata,
//                parameters.w1, parameters.w2, parameters.raised_w);
//
//        add(panelXmippBrowser, BorderLayout.CENTER);
//        pack();
    }

    @Override
    protected JPanelXmippBrowser createPanel(Parameters parameters) {
        return new JPanelXmippBandPassFilter(parameters.files[0],
                parameters.w1, parameters.w2, parameters.raised_w);
    }

    @Override
    protected boolean sendSelectedFiles() {
        JPanelXmippBandPassFilter panel = (JPanelXmippBandPassFilter) panelXmippBrowser;

        return send(new Object[]{panel.getW1(), panel.getW2(), panel.getRaisedW()}, true);
    }
}
