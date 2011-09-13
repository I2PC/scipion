
import browser.COMMAND_PARAMETERS;
import browser.DEBUG;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import wizards.JFrameXmippBandPassFilter;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippBandPassFilterWizard {

    final static int N_ARGS = 5;
    final static int METADATA = 0;
    final static int PORT = 1;
    final static int W1 = 2;
    final static int W2 = 3;
    final static int RAISED_W = 4;

    public static void main(String args[]) {
        String parameters[] = processArgs(args);

        String metadata = parameters[METADATA];
        int port = Integer.parseInt(parameters[PORT]);
        double w1 = Double.parseDouble(parameters[W1]);
        double w2 = Double.parseDouble(parameters[W2]);
        double raised_w = Double.parseDouble(parameters[RAISED_W]);

        JFrameXmippBandPassFilter frameBrowser = new JFrameXmippBandPassFilter(
                metadata, port, w1, w2, raised_w);
        frameBrowser.setVisible(true);
    }

    static String[] processArgs(String args[]) {
        String parameters[] = new String[N_ARGS];
        Options options = new Options();

        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_FILTER_W1, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_FILTER_W2, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_FILTER_RAISED_W, true, "");

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

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_FILTER_W2)) {
                parameters[W2] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_FILTER_W2);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_FILTER_RAISED_W)) {
                parameters[RAISED_W] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_FILTER_RAISED_W);
            }
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }

        return parameters;
    }
}
