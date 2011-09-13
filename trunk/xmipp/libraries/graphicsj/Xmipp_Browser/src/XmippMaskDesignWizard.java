
import browser.COMMAND_PARAMETERS;
import browser.DEBUG;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import wizards.JFrameXmippMaskDesign;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippMaskDesignWizard {

    final static int N_ARGS = 3;
    final static int METADATA = 0;
    final static int PORT = 1;
    final static int MASKSFILENAME = 2;

    public static void main(String args[]) {
        String metadata = "/home/juanjo/temp/angles_save.sel";
        int port = 12345;
        String maskfilename = "/home/juanjo/Desktop/mask_filename_from_wizard.xmp";

        JFrameXmippMaskDesign frameBrowser = new JFrameXmippMaskDesign(metadata,
                port, maskfilename);
        frameBrowser.setVisible(true);
    }

    public static void main_(String args[]) {
        String parameters[] = processArgs(args);

        String metadata = parameters[METADATA];
        int port = Integer.parseInt(parameters[PORT]);
        String maskfilename = parameters[MASKSFILENAME];

        JFrameXmippMaskDesign frameBrowser = new JFrameXmippMaskDesign(metadata,
                port, maskfilename);
        frameBrowser.setVisible(true);
    }

    static String[] processArgs(String args[]) {
        String parameters[] = new String[N_ARGS];
        Options options = new Options();

        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_FILTER_BAD_PIXELS_FACTOR, true, "");

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, args);

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE)) {
                parameters[METADATA] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_INPUT_FILE);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT)) {
                parameters[PORT] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_SOCKET_PORT);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_MASKFILENAME)) {
                parameters[MASKSFILENAME] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_MASKFILENAME);
            }
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }

        return parameters;
    }
}
