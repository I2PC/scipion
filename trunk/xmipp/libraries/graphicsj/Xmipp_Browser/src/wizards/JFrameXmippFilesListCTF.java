/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.COMMAND_PARAMETERS;
import browser.LABELS;
import browser.filebrowsers.JDialogXmippFilesList;
import browser.windows.ImagesWindowFactory;
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

//        // Hack: Replaces panel.
//        remove(panelXmippBrowser);
//
//        panelXmippBrowser = new JPanelXmippFileListCTF(directory, expression, downsampling);
//
//        add(panelXmippBrowser, BorderLayout.CENTER);
//        pack();
        setPanel(directory, expression, false, downsampling);
    }

    // Hack: Replaces panel avoiding the super class to add the old one before.
    @Override
    protected void setPanel(final String directory, final String expression, final boolean singleSelection) {
    }

    void setPanel(final String directory, final String expression, final boolean singleSelection, final double downsampling) {
        ImagesWindowFactory.blockGUI(getRootPane(), "Building list...");

        Thread t = new Thread(new Runnable() {

            public void run() {
                panelXmippBrowser = new JPanelXmippFileListCTF(directory, expression, downsampling);

                add(panelXmippBrowser, BorderLayout.CENTER);

                setLocationRelativeTo(null);
                ImagesWindowFactory.releaseGUI(getRootPane());
            }
        });

        t.start();
    }

    @Override
    protected boolean sendSelectedFiles() {
        return send(new Object[]{
                    ((JPanelXmippFileListCTF) panelXmippBrowser).getDownsampling()}, true);
    }
}
