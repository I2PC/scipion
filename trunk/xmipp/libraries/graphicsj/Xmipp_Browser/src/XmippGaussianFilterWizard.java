
import browser.COMMAND_PARAMETERS;
import browser.DEBUG;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import wizards.JDialogXmippGaussianFilter;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippGaussianFilterWizard {

    final static int N_ARGS = 3;
    final static int METADATA = 0;
    final static int PORT = 1;
    final static int W1 = 2;

    public static void main(String args[]) {
        String parameters[] = processArgs(args);

        String metadata = parameters[METADATA];
        int port = Integer.parseInt(parameters[PORT]);
        double w1 = Double.parseDouble(parameters[W1]);

        JDialogXmippGaussianFilter frameBrowser = new JDialogXmippGaussianFilter(
                metadata, port, w1);
        frameBrowser.setVisible(true);
    }

    static String[] processArgs(String args[]) {
        String parameters[] = new String[N_ARGS];
        Options options = new Options();

        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE, true, COMMAND_PARAMETERS.OPTION_INPUT_FILE_DESCRIPTION);
        options.addOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT, true, COMMAND_PARAMETERS.OPTION_SOCKET_PORT_DESCRIPTION);
        options.addOption(COMMAND_PARAMETERS.OPTION_FILTER_W1, true, COMMAND_PARAMETERS.OPTION_FILTER_W1_DESCRIPTION);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, args);

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE)) {
                parameters[METADATA] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_INPUT_FILE);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT)) {
                parameters[PORT] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_SOCKET_PORT);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_FILTER_W1)) {
                parameters[W1] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_FILTER_W1);
            }

        } catch (Exception ex) {
            DEBUG.printException(ex);
        }

        return parameters;
    }
}
