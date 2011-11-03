/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.LABELS;
import browser.commandline.Parameters;
import browser.filebrowsers.JDialogXmippFilesList;
import browser.filebrowsers.JPanelXmippBrowser;
import java.awt.BorderLayout;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameXmippBadPixelsFilter extends JDialogXmippFilesList {

    public JFrameXmippBadPixelsFilter(String metadata, Parameters parameters) {
        super(metadata, parameters);

        setTitle(LABELS.TITLE_WIZARD_BAD_PIXELS_FILTER);
/*
        // Hack: Replaces panel.
        remove(panelXmippBrowser);

        panelXmippBrowser = new JPanelXmippBadPixelsFilter(metadata, parameters.bad_pixels_factor);

        add(panelXmippBrowser, BorderLayout.CENTER);
        pack();*/
    }

    @Override
    protected JPanelXmippBrowser createPanel(Parameters parameters) {
        return new JPanelXmippBadPixelsFilter(parameters.files[0], parameters.bad_pixels_factor);
    }

    @Override
    protected boolean sendSelectedFiles() {
        JPanelXmippBadPixelsFilter panel = (JPanelXmippBadPixelsFilter) panelXmippBrowser;

        return send(new Object[]{panel.getFactor()}, true);
    }
}
