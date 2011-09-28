
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
public class Stitching {

    final static String OPTION_PROPERTIES_FILE = "p";
    final static String OPTION_INPUT_FILE = "i";
    final static String OPTION_OUTPUT_FILE = "o";
    String FILES[];
    String PROPERTIES_FILE;
    String OUTPUT_FILE;

    public static void main(String args[]) {
        Stitching stitching = new Stitching(args);
    }

    public Stitching(String args[]) {
        Stitcher stitcher;

        processArgs(args);

        stitcher = new Stitcher(FILES, PROPERTIES_FILE, OUTPUT_FILE);
        stitcher.run();
    }

    final void processArgs(String args[]) {
        Options options = new Options();

        options.addOption(OPTION_INPUT_FILE, true, "");
        options.addOption(OPTION_PROPERTIES_FILE, true, "");
        options.addOption(OPTION_OUTPUT_FILE, true, "");

        // It should be able to handle multiple files.
        options.getOption(OPTION_INPUT_FILE).setArgs(Integer.MAX_VALUE);
//        options.getOption(OPTION_PROPERTIES_FILE).setOptionalArg(true);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, args);

            if (cmdLine.hasOption(OPTION_INPUT_FILE)) {
                FILES = cmdLine.getOptionValues(OPTION_INPUT_FILE);
            }

            if (cmdLine.hasOption(OPTION_PROPERTIES_FILE)) {
                PROPERTIES_FILE = cmdLine.getOptionValue(OPTION_PROPERTIES_FILE);
            }

            if (cmdLine.hasOption(OPTION_OUTPUT_FILE)) {
                OUTPUT_FILE = cmdLine.getOptionValue(OPTION_OUTPUT_FILE);
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }
}
