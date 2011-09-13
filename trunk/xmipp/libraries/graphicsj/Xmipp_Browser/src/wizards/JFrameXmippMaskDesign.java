/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.LABELS;
import browser.filebrowsers.JDialogXmippFilesList;
import browser.imageitems.ImageConverter;
import ij.IJ;
import ij.ImagePlus;
import java.awt.BorderLayout;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameXmippMaskDesign extends JDialogXmippFilesList {

    String maskfilename;

    public JFrameXmippMaskDesign(String metadata, int port, String maskfilename) {
        super(metadata, port);

        this.maskfilename = maskfilename;

        setTitle(LABELS.TITLE_WIZARD_MASK_DESIGN);

        // Hack: Replaces panel.
        remove(panelXmippBrowser);

        panelXmippBrowser = new JPanelXmippMaskDesign(metadata);

        add(panelXmippBrowser, BorderLayout.CENTER);

        pack();
    }

    @Override
    protected boolean sendSelectedFiles() {
        ImagePlus imp = IJ.getImage();
        IJ.run(imp, "32-bit", "");

        if (ImageConverter.saveImage(imp, maskfilename)) {
            IJ.showMessage("Mask saved!: " + maskfilename);
        }

        return send(null, true);
    }
}
