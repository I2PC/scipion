
import browser.COMMAND_PARAMETERS;
import browser.filebrowsers.JFrameXmippFilesList;
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
public class XmippFileList extends XmippBrowser {

    int PORT;

    public XmippFileList(int PORT) {
        super();

        this.PORT = PORT;
    }

    @Override
    void runBrowser(String directory, String expression, boolean singleSelection) {
        JFrameXmippFilesList frameBrowser = new JFrameXmippFilesList(
                directory, PORT, expression, singleSelection);
        frameBrowser.setVisible(true);
    }

    @Override
    void processArgs(String args) {
        super.processArgs(args);

        String argsList[] = args.split(" ");
        Options options = new Options();

        options.addOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT, true, COMMAND_PARAMETERS.OPTION_SOCKET_PORT_DESCRIPTION);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_SOCKET_PORT)) {
                PORT = Integer.valueOf(cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_SOCKET_PORT));
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }
}
