
import browser.COMMAND_PARAMETERS;
import browser.DEBUG;
import java.text.SimpleDateFormat;
import java.util.Date;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import wizards.JFrameXmippFilesListCTF;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippFileListCTFWizard {

    final static int N_ARGS = 4;
    final static int DIR = 0;
    final static int PORT = 1;
//    final static int SINGLESEL = 2;
//    final static int SELTYPE = 3;
    final static int FILTER = 2;
    final static int DOWNSAMPLING = 3;

    public static void main(String args[]) {
        long before = System.currentTimeMillis();
        String parameters[] = processArgs(args);

        String dir = parameters[DIR];
        int port = Integer.parseInt(parameters[PORT]);
        //boolean singleSel = Boolean.parseBoolean(parameters[SINGLESEL]);
        //String selType = parameters[SELTYPE];
        String expression = parameters[FILTER];
        double downsampling = Double.parseDouble(parameters[DOWNSAMPLING]);

        JFrameXmippFilesListCTF frameBrowser = new JFrameXmippFilesListCTF(
                dir, port, expression, downsampling);
        frameBrowser.setVisible(true);

        long after = System.currentTimeMillis();
        long elapsed = after - before;

        SimpleDateFormat dateFormatter = new SimpleDateFormat("'Time to start:' mm:ss:S");
        System.out.println(dateFormatter.format(new Date(elapsed)));
    }

    static String[] processArgs(String args[]) {
        String parameters[] = new String[N_ARGS];
        Options options = new Options();

        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_DIR, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_FILTER, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_FILTER_DOWNSAMPLING, true, "");

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, args);

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_DIR)) {
                parameters[DIR] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_INPUT_DIR);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT)) {
                parameters[PORT] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_SOCKET_PORT);
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

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_FILTER_DOWNSAMPLING)) {
                parameters[DOWNSAMPLING] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_FILTER_DOWNSAMPLING);
            }
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }

        return parameters;
    }
}
