/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.commandline;

import browser.DEBUG;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

/**
 *
 * @author Juanjo Vega
 */
public class Parameters {

    public String files[];
    public int port;
    public String filter;
    public String mode;
    public boolean poll;
    public int zoom = 100;
    public boolean renderImages;
    public int rows = -1, columns = -1;
    public String vectorsFilename;
    public String classesFilename;
    public String dataFilename;

    public Parameters() {
    }

    public Parameters(String args[]) {
        processArgs(args);
    }

    void processArgs(String args[]) {
        Options options = new Options();

        options.addOption(COMMAND_PARAMETERS.FILE, true, "");
        // It should be able to handle multiple files.
        options.getOption(COMMAND_PARAMETERS.FILE).setOptionalArg(true);
        options.getOption(COMMAND_PARAMETERS.FILE).setArgs(Integer.MAX_VALUE);

        options.addOption(COMMAND_PARAMETERS.OPENING_MODE, true, "");
        options.addOption(COMMAND_PARAMETERS.POLL, false, "");
        options.addOption(COMMAND_PARAMETERS.ZOOM, true, "");
        options.addOption(COMMAND_PARAMETERS.RENDER_IMAGES, false, "");
        options.addOption(COMMAND_PARAMETERS.TABLE_ROWS, true, "");
        options.addOption(COMMAND_PARAMETERS.TABLE_COLUMNS, true, "");

        options.addOption(COMMAND_PARAMETERS.INPUT_VECTORSFILE, true, "");
        options.addOption(COMMAND_PARAMETERS.INPUT_CLASSESFILE, true, "");
        options.addOption(COMMAND_PARAMETERS.INPUT_DATAFILE, true, "");

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, args);

            if (cmdLine.hasOption(COMMAND_PARAMETERS.FILE)) {
                files = cmdLine.getOptionValues(COMMAND_PARAMETERS.FILE);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPENING_MODE)) {
                mode = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPENING_MODE);
            } else {
                mode = COMMAND_PARAMETERS.OPENING_MODE_DEFAULT;
            }

            poll = cmdLine.hasOption(COMMAND_PARAMETERS.POLL);

            if (cmdLine.hasOption(COMMAND_PARAMETERS.ZOOM)) {
                zoom = Integer.parseInt(cmdLine.getOptionValue(COMMAND_PARAMETERS.ZOOM));
            }

            renderImages = cmdLine.hasOption(COMMAND_PARAMETERS.RENDER_IMAGES);

            if (cmdLine.hasOption(COMMAND_PARAMETERS.TABLE_ROWS)) {
                rows = Integer.parseInt(cmdLine.getOptionValue(COMMAND_PARAMETERS.TABLE_ROWS));
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.TABLE_COLUMNS)) {
                columns = Integer.parseInt(cmdLine.getOptionValue(COMMAND_PARAMETERS.TABLE_COLUMNS));
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.INPUT_VECTORSFILE)) {
                vectorsFilename = cmdLine.getOptionValue(COMMAND_PARAMETERS.INPUT_VECTORSFILE);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.INPUT_CLASSESFILE)) {
                classesFilename = cmdLine.getOptionValue(COMMAND_PARAMETERS.INPUT_CLASSESFILE);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.INPUT_DATAFILE)) {
                dataFilename = cmdLine.getOptionValue(COMMAND_PARAMETERS.INPUT_DATAFILE);
            }
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }
    }
}
