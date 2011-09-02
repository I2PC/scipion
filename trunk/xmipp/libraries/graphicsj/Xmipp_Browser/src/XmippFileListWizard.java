
import browser.COMMAND_PARAMETERS;
import browser.DEBUG;
import browser.filebrowsers.JDialogXmippFilesList;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippFileListWizard {

    final static int N_ARGS = 5;
    final static int DIR = 0;
    final static int PORT = 1;
    final static int SINGLESEL = 2;
    final static int SELTYPE = 3;
    final static int FILTER = 4;

    public static void main(String args[]) {
        String parameters[] = processArgs(args);

        String dir = parameters[DIR];
        int port = Integer.parseInt(parameters[PORT]);
        boolean singleSel = Boolean.parseBoolean(parameters[SINGLESEL]);
        String selType = parameters[SELTYPE];
        String filterExp = parameters[FILTER];

        JDialogXmippFilesList frameBrowser = new JDialogXmippFilesList(
                dir, port, singleSel, selType, filterExp);
        frameBrowser.setVisible(true);
    }

    static String[] processArgs(String args[]) {
        String parameters[] = new String[N_ARGS];
        Options options = new Options();

        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_DIR, true, COMMAND_PARAMETERS.OPTION_INPUT_DIR_DESCRIPTION);
        options.addOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT, true, COMMAND_PARAMETERS.OPTION_SOCKET_PORT_DESCRIPTION);
        options.addOption(COMMAND_PARAMETERS.OPTION_SINGLE_SELECTION, false, COMMAND_PARAMETERS.OPTION_SINGLE_SELECTION_DESCRIPTION);
        options.addOption(COMMAND_PARAMETERS.OPTION_SELECTION_TYPE, true, COMMAND_PARAMETERS.OPTION_SELECTION_TYPE_DESCRIPTION);
        options.addOption(COMMAND_PARAMETERS.OPTION_FILTER, true, COMMAND_PARAMETERS.OPTION_FILTER_DESCRIPTION);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, args);

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_DIR)) {
                parameters[DIR] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_INPUT_DIR);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT)) {
                parameters[PORT] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_SOCKET_PORT);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_SINGLE_SELECTION)) {
                parameters[SINGLESEL] = "true";
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_SINGLE_SELECTION)) {
                parameters[SELTYPE] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_SELECTION_TYPE);
            } else {
                parameters[SELTYPE] = COMMAND_PARAMETERS.SELECTION_TYPE_ANY;
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_FILTER)) {
                String filters[] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_FILTER).split(COMMAND_PARAMETERS.FILTERS_SEPARATOR);

                String filter = "";
                for (int i = 0; i < filters.length; i++) {
                    filter += filters[i];
                    if (i < filters.length - 1) {
                        filter += " ";
                    }
                }
                parameters[FILTER] = filter;
            }

        } catch (Exception ex) {
            DEBUG.printException(ex);
        }

        return parameters;
    }
}
