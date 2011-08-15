
import ij.IJ;
import ij.ImageJ;
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
public class Stitching implements PlugIn {

    private final static String OPTION_INPUT_FILE = "i";
    private final static String OPTION_INPUT_FILE_DESCRIPTION = "input file(s)";
    String FILES[];

    public static void main(String args[]) {
        new ImageJ();

        String basedir = "/home/juanjo/Desktop/";
        String files[] = {
            basedir + "Clipboard2.tif",
            basedir + "Clipboard1.tif"
        };

        Stitcher stitcher = new Stitcher(files);
        stitcher.run();
    }

    public void run(String string) {
        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "string" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            processArgs(Macro.getOptions().split(" "));

            Stitcher stitcher = new Stitcher(FILES);
            stitcher.run();
        } else {    // From menu.
            IJ.error("No menu run implementation yet.");
        }
    }

    private void processArgs(String args[]) {
        Options options = new Options();

        options.addOption(OPTION_INPUT_FILE, true, OPTION_INPUT_FILE_DESCRIPTION);

        // It should be able to handle multiple files.
        options.getOption(OPTION_INPUT_FILE).setOptionalArg(true);
        options.getOption(OPTION_INPUT_FILE).setArgs(Integer.MAX_VALUE);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, args);

            if (cmdLine.hasOption(OPTION_INPUT_FILE)) {
                FILES = cmdLine.getOptionValues(OPTION_INPUT_FILE);
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }
}
