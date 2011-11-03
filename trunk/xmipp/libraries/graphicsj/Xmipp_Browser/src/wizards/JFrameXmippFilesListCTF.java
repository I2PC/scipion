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
public class JFrameXmippFilesListCTF extends JDialogXmippFilesList {

    public JFrameXmippFilesListCTF(String directory, Parameters parameters) {
        super(directory, parameters);

        setTitle(LABELS.TITLE_WIZARD_PSD);

//        // Hack: Replaces panel.
//        remove(panelXmippBrowser);
//
//        panelXmippBrowser = new JPanelXmippFileListCTF(directory, expression, downsampling);
//
//        add(panelXmippBrowser, BorderLayout.CENTER);
//        pack();
//        setPanel(parameters);
    }

    // Hack: Replaces panel avoiding the super class to add the old one before.
    @Override
    protected JPanelXmippBrowser createPanel(Parameters parameters) {
//    protected JPanelXmippBrowser setPanel(final String directory, final Parameters parameters) {
        panelXmippBrowser = new JPanelXmippFileListCTF(parameters.directory,
                parameters.filter, parameters.downsampling);

        return panelXmippBrowser;
    }

    @Override
    protected boolean sendSelectedFiles() {
        return send(new Object[]{
                    ((JPanelXmippFileListCTF) panelXmippBrowser).getDownsampling()}, true);
    }
}
