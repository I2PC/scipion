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
public class JDialogXmippFilesListMeasureFreqsBrowser extends JDialogXmippFilesListCTF {

    public JDialogXmippFilesListMeasureFreqsBrowser(String directory, int port, String expression) {
        this(directory, port, expression, 1.0);
    }

    public JDialogXmippFilesListMeasureFreqsBrowser(String directory, int port, String expression, double downsampling) {
        super(directory, port, expression);

        setTitle(LABELS.TITLE_XMIPP_WIZARD_MEASURE_FREQS);

        jbOk.setText(LABELS.BUTTON_CLOSE);
        jpButtons.remove(jbCancel);

        // Hack: Replaces panel.
        remove(panelXmippBrowser);

        panelXmippBrowser = new JPanelXmippFileListCTF(directory, expression, downsampling);
        ((JPanelXmippFileListCTF) panelXmippBrowser).setDownsamplingEnabled(false);

        add(panelXmippBrowser, BorderLayout.CENTER);
        pack();
    }

    @Override
    protected void button1Clicked() {
        cancel();
    }

    @Override
    void cancel() {
        dispose();
    }
}
