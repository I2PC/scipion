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
public class JFrameXmippGaussianFilter extends JDialogXmippFilesList {

    public JFrameXmippGaussianFilter(String metadata, Parameters parameters) {
        super("", parameters);

        setTitle(LABELS.TITLE_WIZARD_GAUSSIAN_FILTER);
//
//        // Hack: Replaces panel.
//        remove(panelXmippBrowser);
//
//        panelXmippBrowser = new JPanelXmippGaussianFilter(metadata, parameters.w1);
//
//        add(panelXmippBrowser, BorderLayout.CENTER);
//        pack();
    }

    @Override
    protected JPanelXmippBrowser createPanel(Parameters parameters) {
        return new JPanelXmippGaussianFilter(parameters.files[0], parameters.w1);
    }

    @Override
    protected boolean sendSelectedFiles() {
        JPanelXmippGaussianFilter panel = (JPanelXmippGaussianFilter) panelXmippBrowser;

        return send(new Object[]{panel.getW1()}, true);
    }
}
