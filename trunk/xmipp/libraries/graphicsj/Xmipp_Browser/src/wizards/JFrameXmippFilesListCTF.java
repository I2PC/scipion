/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.COMMAND_PARAMETERS;
import browser.LABELS;
import browser.filebrowsers.JDialogXmippFilesList;
import java.awt.BorderLayout;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameXmippFilesListCTF extends JDialogXmippFilesList {

    public JFrameXmippFilesListCTF(String directory, int port, String expression) {
        this(directory, port, expression, 1.0);
    }

    public JFrameXmippFilesListCTF(String directory, int port, String expression, double downsampling) {
        super(directory, port, false, COMMAND_PARAMETERS.SELECTION_TYPE_ANY, expression);

        setTitle(LABELS.TITLE_WIZARD_PSD);

        // Hack: Replaces panel.
        remove(panelXmippBrowser);

        panelXmippBrowser = new JPanelXmippFileListCTF(directory, expression, downsampling);

        add(panelXmippBrowser, BorderLayout.CENTER);
        pack();
    }

    @Override
    protected boolean sendSelectedFiles() {
        return send(new Object[]{
                    ((JPanelXmippFileListCTF) panelXmippBrowser).getDownsampling()}, true);
    }
}
