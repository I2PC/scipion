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
public class JFrameXmippBandPassFilter extends JDialogXmippFilesList {

    public JFrameXmippBandPassFilter(String metadata, int port) {
        this(metadata, port, 1.0, 1.0, 1.0);
    }

    public JFrameXmippBandPassFilter(String metadata, int port, double w1, double w2, double raised_w) {
        super(metadata, port);

        setTitle(LABELS.TITLE_WIZARD_BAND_PASS_FILTER);

        // Hack: Replaces panel.
        remove(panelXmippBrowser);

        panelXmippBrowser = new JPanelXmippBandPassFilter(metadata, w1, w2, raised_w);

        add(panelXmippBrowser, BorderLayout.CENTER);
        pack();
    }

    @Override
    protected boolean sendSelectedFiles() {
        JPanelXmippBandPassFilter panel = (JPanelXmippBandPassFilter) panelXmippBrowser;

        return send(new Object[]{panel.getW1(), panel.getW2(), panel.getRaisedW()}, true);
    }
}
