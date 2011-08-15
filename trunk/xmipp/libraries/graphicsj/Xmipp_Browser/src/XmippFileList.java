
import browser.COMMAND_PARAMETERS;
import browser.filebrowsers.JDialogXmippFilesList;
import ij.IJ;
import ij.Macro;
import ij.plugin.PlugIn;
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
public class XmippFileList implements PlugIn {

    // Browser
    private String DIR;
    private boolean SINGLE_SELECTION = false;
    private String FILTER = "";
    int PORT;

    @Override
    public void run(String string) {
        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "string" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            processArgs(Macro.getOptions().trim());
        }

        if (DIR == null) {
            DIR = System.getProperty("user.dir");
        }

        runBrowser(DIR, PORT, FILTER, SINGLE_SELECTION);
    }

    void runBrowser(String directory, int port, String expression, boolean singleSelection) {
//        IJ.getInstance().setExtendedState(Frame.ICONIFIED);
//        IJ.getInstance().setVisible(false);

        JDialogXmippFilesList frameBrowser = new JDialogXmippFilesList(directory, port, expression, singleSelection);
        frameBrowser.setVisible(true);

//        IJ.getInstance().setVisible(true);
    }

    void processArgs(String args) {
        String argsList[] = args.split(" ");
        Options options = new Options();

        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_DIR, true, COMMAND_PARAMETERS.OPTION_INPUT_DIR_DESCRIPTION);
        options.addOption(COMMAND_PARAMETERS.OPTION_FILTER, true, COMMAND_PARAMETERS.OPTION_FILTER);
        options.addOption(COMMAND_PARAMETERS.OPTION_SINGLE_SELECTION, false, COMMAND_PARAMETERS.OPTION_SINGLE_SELECTION);

        options.addOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT, true, COMMAND_PARAMETERS.OPTION_SOCKET_PORT_DESCRIPTION);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_DIR)) {
                DIR = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_INPUT_DIR);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_FILTER)) {
                String filters[] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_FILTER).split(COMMAND_PARAMETERS.FILTERS_SEPARATOR);

                FILTER = "";
                for (int i = 0; i < filters.length; i++) {
                    FILTER += filters[i];
                    if (i < filters.length - 1) {
                        FILTER += " ";
                    }
                }
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_SINGLE_SELECTION)) {
                SINGLE_SELECTION = true;
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT)) {
                PORT = Integer.valueOf(cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_SOCKET_PORT));
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }
}
