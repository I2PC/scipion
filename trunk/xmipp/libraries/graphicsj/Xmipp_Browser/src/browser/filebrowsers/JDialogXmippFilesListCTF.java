/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.filebrowsers;

import browser.LABELS;
import java.awt.BorderLayout;

/**
 *
 * @author Juanjo Vega
 */
public class JDialogXmippFilesListCTF extends JDialogXmippFilesList {

    /** Creates new form JDialogXmippFilesList */
    public JDialogXmippFilesListCTF(String directory, int port) {
        this(directory, port, false);
    }

    public JDialogXmippFilesListCTF(String directory, int port, String expression) {
        this(directory, port, expression, false);
    }

    public JDialogXmippFilesListCTF(String directory, int port, boolean singleSelection) {
        this(directory, port, "", singleSelection);
    }

    public JDialogXmippFilesListCTF(String directory, int port, String expression, boolean singleSelection) {
        super(directory, port, expression, singleSelection);

        setTitle(LABELS.TITLE_XMIPP_FILE_SELECTOR_CTF);

        // Hack: Replaces panel.
        remove(panelXmippBrowser);

        panelXmippBrowser = new JPanelXmippFileListCTF(directory, expression);
        panelXmippBrowser.setSingleSelection(singleSelection);

        add(panelXmippBrowser, BorderLayout.CENTER);
    }
}
