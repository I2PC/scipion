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
public class JFrameXmippBadPixelsFilter extends JDialogXmippFilesList {

    public JFrameXmippBadPixelsFilter(String metadata, int port) {
        this(metadata, port, 1.0);
    }

    public JFrameXmippBadPixelsFilter(String metadata, int port, double factor) {
        super(metadata, port);

        setTitle(LABELS.TITLE_WIZARD_BAD_PIXELS_FILTER);

        // Hack: Replaces panel.
        remove(panelXmippBrowser);

        panelXmippBrowser = new JPanelXmippBadPixelsFilter(metadata, factor);

        add(panelXmippBrowser, BorderLayout.CENTER);
        pack();
    }

    @Override
    protected boolean sendSelectedFiles() {
        JPanelXmippBadPixelsFilter panel = (JPanelXmippBadPixelsFilter) panelXmippBrowser;

        return send(new Object[]{panel.getFactor()}, true);
    }
}
