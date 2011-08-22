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
//    public JDialogXmippFilesListCTF(String directory, int port) {
//        this(directory, port, false);
//    }
//
//    public JDialogXmippFilesListCTF(String directory, int port, String expression) {
//        this(directory, port, expression, false);
//    }
//    public JDialogXmippFilesListCTF(String directory, int port, boolean singleSelection) {
//        this(directory, port, "", singleSelection);
//    }
    public JDialogXmippFilesListCTF(String directory, int port, String expression) {
        this(directory, port, expression, 1.0);
    }

    public JDialogXmippFilesListCTF(String directory, int port, String expression, double downsampling) {
        super(directory, port, expression);

        setTitle(LABELS.TITLE_XMIPP_WIZARD_PSD);

        // Hack: Replaces panel.
        remove(panelXmippBrowser);

        panelXmippBrowser = new JPanelXmippFileListCTF(directory, expression, downsampling);

        add(panelXmippBrowser, BorderLayout.CENTER);
        pack();
    }

    @Override
    protected boolean sendSelectedFiles() {
        return send(new Object[]{
                    ((JPanelXmippFileListCTF) panelXmippBrowser).getDownsampling()});
    }
}
