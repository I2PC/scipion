

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import browser.COMMAND_PARAMETERS;
import browser.DEBUG;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import ij.Macro;
import ij.plugin.PlugIn;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

/**
 *
 * @author Juanjo Vega
 */
public class XmippRotSpectraViewer implements PlugIn {

    String VECTORS, CLASSES, DATA;

    public void run(String args) {
        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "args" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            processArgs(Macro.getOptions().trim());

            if (VECTORS != null && CLASSES != null && DATA != null) {
                openFile(VECTORS, CLASSES, DATA);
            }
        }
    }

    public static void openFile(String vectorsFile, String classesFile, String dataFile) {
        ImagesWindowFactory.openRotSpectrasWindow(vectorsFile, classesFile, dataFile);
    }

    void processArgs(String args) {
        String argsList[] = args.split(" ");
        Options options = new Options();

        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_VECTORSFILE, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_CLASSESFILE, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_DATAFILE, true, "");

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            // Input.
            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_VECTORSFILE)) {
                VECTORS = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_INPUT_VECTORSFILE);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_CLASSESFILE)) {
                CLASSES = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_INPUT_CLASSESFILE);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_DATAFILE)) {
                DATA = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_INPUT_DATAFILE);
            }

        } catch (Exception ex) {
            DEBUG.printException(ex);
            throw new RuntimeException(ex);
        }
    }
}
