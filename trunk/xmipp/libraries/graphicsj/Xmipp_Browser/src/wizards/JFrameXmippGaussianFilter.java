/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.LABELS;
import browser.filebrowsers.JDialogXmippFilesList;
import java.awt.BorderLayout;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameXmippGaussianFilter extends JDialogXmippFilesList {

    public JFrameXmippGaussianFilter(String metadata, int port, double w1) {
        super("", port);

        setTitle(LABELS.TITLE_WIZARD_GAUSSIAN_FILTER);

        // Hack: Replaces panel.
        remove(panelXmippBrowser);

        panelXmippBrowser = new JPanelXmippGaussianFilter(metadata, w1);

        add(panelXmippBrowser, BorderLayout.CENTER);
        pack();
    }

    @Override
    protected boolean sendSelectedFiles() {
        JPanelXmippGaussianFilter panel = (JPanelXmippGaussianFilter) panelXmippBrowser;

        return send(new Object[]{panel.getW1()}, true);
    }
}
