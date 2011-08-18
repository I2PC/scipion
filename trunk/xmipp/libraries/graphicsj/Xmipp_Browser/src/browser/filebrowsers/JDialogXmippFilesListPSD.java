/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.filebrowsers;

import java.awt.BorderLayout;

/**
 *
 * @author Juanjo Vega
 */
public class JDialogXmippFilesListPSD extends JDialogXmippFilesList {

    /** Creates new form JDialogXmippFilesList */
    public JDialogXmippFilesListPSD(String directory, int port) {
        this(directory, port, false);
    }

    public JDialogXmippFilesListPSD(String directory, int port, String expression) {
        this(directory, port, expression, false);
    }

    public JDialogXmippFilesListPSD(String directory, int port, boolean singleSelection) {
        this(directory, port, "", singleSelection);
    }

    public JDialogXmippFilesListPSD(String directory, int port, String expression, boolean singleSelection) {
        super(directory, port, expression, singleSelection);

        // Hack: Replaces panel.
        remove(panelXmippBrowser);

        panelXmippBrowser = new JPanelXmippFileListPSD(directory, expression);
        panelXmippBrowser.setSingleSelection(singleSelection);

        add(panelXmippBrowser, BorderLayout.CENTER);
    }

    @Override
    protected boolean sendSelectedFiles() {
        return send(new Object[]{
                    ((JPanelXmippFileListPSD) panelXmippBrowser).getDownsampling()});
    }
}
