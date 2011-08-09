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
public class JFrameXmippFilesListPSD extends JFrameXmippFilesList {

    /** Creates new form JFrameXmippFilesList */
    public JFrameXmippFilesListPSD(String directory, int port) {
        this(directory, port, false);
    }

    public JFrameXmippFilesListPSD(String directory, int port, String expression) {
        this(directory, port, expression, false);
    }

    public JFrameXmippFilesListPSD(String directory, int port, boolean singleSelection) {
        this(directory, port, "", singleSelection);
    }

    public JFrameXmippFilesListPSD(String directory, int port, String expression, boolean singleSelection) {
        super(directory, port, expression, singleSelection);

        setTitle(LABELS.TITLE_XMIPP_FILE_SELECTOR_PSD);

        // Hack: Replaces panel.
        remove(panelXmippBrowser);

        panelXmippBrowser = new JPanelXmippFileListPSD(directory, expression);
        panelXmippBrowser.setSingleSelection(singleSelection);

        add(panelXmippBrowser, BorderLayout.CENTER);
    }
}
