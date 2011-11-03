

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import browser.commandline.Parameters;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import ij.Macro;
import ij.plugin.PlugIn;

/**
 *
 * @author Juanjo Vega
 */
public class XmippRotSpectraViewer implements PlugIn {

    public void run(String args) {
        Parameters parameters;

        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "args" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            parameters = new Parameters(Macro.getOptions().split(" "));

            if (parameters.vectorsFilename != null
                    && parameters.classesFilename != null
                    && parameters.dataFilename != null) {
                openFile(parameters.vectorsFilename,
                        parameters.classesFilename,
                        parameters.dataFilename);
            }
        }
    }

    public static void openFile(String vectorsFile, String classesFile, String dataFile) {
        ImagesWindowFactory.openRotSpectrasWindow(vectorsFile, classesFile, dataFile);
    }
}
