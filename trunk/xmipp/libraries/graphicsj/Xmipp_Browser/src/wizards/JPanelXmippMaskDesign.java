/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.imageitems.listitems.AbstractImageItem;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Toolbar;

/**
 *
 * @author Juanjo Vega
 */
public class JPanelXmippMaskDesign extends JPanelXmippFilterMetadata {

    String maskfilename;

    public JPanelXmippMaskDesign(String metadata, String maskfilename) {
        super(metadata);

        this.maskfilename = maskfilename;
        jpPreview.remove(jlFilter);
    }

    @Override
    ImagePlus getFilteredPreview(AbstractImageItem item) throws Exception {
        return getPreview(item);
    }

    @Override
    protected void openSelectedFile() {
        int currentTool = Toolbar.getToolId();

        super.openSelectedFile();

        IJ.setTool(currentTool);

        IJ.run("Masks Tool Bar", maskfilename);
    }
}
