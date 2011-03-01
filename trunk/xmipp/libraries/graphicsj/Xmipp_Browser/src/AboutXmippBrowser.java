
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
public class AboutXmippBrowser implements PlugIn {

    public void run(String string) {
        IJ.showMessage("About Xmipp Browser plugin...",
                "Juanjo Vega"
                + "\nBiocomputing Unit."
                + "\nNational Center for Biotechnology (CNB/CSIC)."
                + "\nMadrid. Jan 2010.");
    }
}
