/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.LABELS;
import java.awt.BorderLayout;

/**
 *
 * @author Juanjo Vega
 */
public class JDialogXmippFilesListMeasureFreqsBrowser extends JDialogXmippFilesListCTF {

    public JDialogXmippFilesListMeasureFreqsBrowser(String directory, int port, String expression) {
        this(directory, port, expression, 1.0);
    }

    public JDialogXmippFilesListMeasureFreqsBrowser(String directory, int port, String expression, double downsampling) {
        super(directory, port, expression);

        setTitle(LABELS.TITLE_WIZARD_MEASURE_FREQS);

        jpButtons.remove(jbOk);
        jbCancel.setText(LABELS.BUTTON_CLOSE);

        // Hack: Replaces panel.
        remove(panelXmippBrowser);

        panelXmippBrowser = new JPanelXmippFileListCTF(directory, expression, downsampling);
        ((JPanelXmippFileListCTF) panelXmippBrowser).setDownsamplingEnabled(false);

        add(panelXmippBrowser, BorderLayout.CENTER);
        pack();
    }

    @Override
    protected void cancel() {
        dispose();
    }
}
