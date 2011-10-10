
import ij.IJ;
import java.text.SimpleDateFormat;
import java.util.Date;
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

    final static String OPTION_INPUT_FILE = "i";
    final static String OPTION_PROPERTIES_FILE = "p";
    final static String OPTION_OUTPUT_FILE = "o";
    final static String OPTION_GENERATESTACK = "stack";
    final static String OPTION_X = "x";
    final static String OPTION_Y = "y";
    String FILES[];
    String PROPERTIES_FILE;
    String OUTPUT_FILE;
    String STACK_FILE = null;
    int X = 0, Y = 0;

    public static void main(String args[]) {
        String basedir = "/gpfs/fs1/home/bioinfo/jvega/";
        //String basedir = "/home/juanjo/Desktop/";
        String files[] = {
            basedir + "left.tif",
            //basedir + "left_2.tif",
            basedir + "right.tif"
        };
        String properties = "stitching.properties";
        String output = basedir + "stitching.xmp";
        String stackfile = basedir + "stitching.stk";

        String arguments[] = new String[]{
            "-" + OPTION_INPUT_FILE, files[0], files[1],
            "-" + OPTION_PROPERTIES_FILE, properties,
            "-" + OPTION_OUTPUT_FILE, output,
            "-" + OPTION_GENERATESTACK, stackfile,
            //"-" + OPTION_X, "100",
            //"-" + OPTION_Y, "-50"
        };

        Stitching stitching = new Stitching(arguments);
    }

    public Stitching(String args[]) {
        processArgs(args);

        // Check parameters.
        if (FILES.length > 1) {
            if (PROPERTIES_FILE != null) {
                if (OUTPUT_FILE != null) {
                    long before = System.currentTimeMillis();

                    Stitcher stitcher = new Stitcher(FILES, PROPERTIES_FILE, OUTPUT_FILE, X, Y, STACK_FILE);
                    stitcher.run();

                    long elapsed = System.currentTimeMillis() - before;
                    SimpleDateFormat dateFormatter = new SimpleDateFormat("'Total time:' mm:ss:S");

                    System.out.println(dateFormatter.format(new Date(elapsed)));
                } else {
                    IJ.error("ERROR: Please, provide a filename for output.");
                }
            } else {
                IJ.error("ERROR: Please, provide a parameters file.");
            }
        } else {
            IJ.error("ERROR: Please, provide at least two files to stitch.");
        }
    }

    final void processArgs(String args[]) {
        Options options = new Options();

        options.addOption(OPTION_INPUT_FILE, true, "");
        options.addOption(OPTION_PROPERTIES_FILE, true, "");
        options.addOption(OPTION_OUTPUT_FILE, true, "");
        options.addOption(OPTION_GENERATESTACK, true, "");
        options.addOption(OPTION_X, true, "");
        options.addOption(OPTION_Y, true, "");

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

            if (cmdLine.hasOption(OPTION_GENERATESTACK)) {
                STACK_FILE = cmdLine.getOptionValue(OPTION_GENERATESTACK);
            }

            if (cmdLine.hasOption(OPTION_X)) {
                X = Integer.valueOf(cmdLine.getOptionValue(OPTION_X));
            }

            if (cmdLine.hasOption(OPTION_Y)) {
                Y = Integer.valueOf(cmdLine.getOptionValue(OPTION_Y));
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }
}
