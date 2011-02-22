
import ij.IJ;
import ij.plugin.PlugIn;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class About_Xmipp_IO implements PlugIn {

    public void run(String string) {
        IJ.showMessage("About Xmipp IO plugin...",
                "Juanjo Vega"
                + "\nBiocomputing Unit."
                + "\nNational Center for Biotechnology (CNB/CSIC)."
                + "\nMadrid. Feb 2011.");
    }
}
