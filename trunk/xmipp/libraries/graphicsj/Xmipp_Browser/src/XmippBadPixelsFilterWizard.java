


import browser.COMMAND_PARAMETERS;
import browser.DEBUG;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import wizards.JFrameXmippBadPixelsFilter;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippBadPixelsFilterWizard {

    final static int N_ARGS = 3;
    final static int METADATA = 0;
    final static int PORT = 1;
    final static int FACTOR = 2;

    public static void main(String args[]) {
        String parameters[] = processArgs(args);

        String metadata = parameters[METADATA];
        int port = Integer.parseInt(parameters[PORT]);
        double factor = Double.parseDouble(parameters[FACTOR]);

        JFrameXmippBadPixelsFilter frameBrowser = new JFrameXmippBadPixelsFilter(
                metadata, port, factor);
        frameBrowser.setVisible(true);
    }

    static String[] processArgs(String args[]) {
        String parameters[] = new String[N_ARGS];
        Options options = new Options();

        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_FILTER_BAD_PIXELS_FACTOR, true, COMMAND_PARAMETERS.OPTION_FILTER_BAD_PIXELS_FACTOR);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, args);

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE)) {
                parameters[METADATA] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_INPUT_FILE);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT)) {
                parameters[PORT] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_SOCKET_PORT);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_FILTER_BAD_PIXELS_FACTOR)) {
                parameters[FACTOR] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_FILTER_BAD_PIXELS_FACTOR);
            }

        } catch (Exception ex) {
            DEBUG.printException(ex);
        }

        return parameters;
    }
}
