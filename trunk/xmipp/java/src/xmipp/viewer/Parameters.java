package xmipp.viewer;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;


/**
 *
 * @author Juanjo Vega
 */
public class Parameters {

    public String directory;
    public String files[];
    public int port;
    public String filter;
    public boolean singleSelection;
    public String selectionType = COMMAND_PARAMETERS.SELECTION_TYPE_ANY;
    public String mode;
    public boolean poll;
    public int zoom = 100;
    public boolean renderImages;
    public int rows = -1, columns = -1;
    public double bad_pixels_factor;
    public double w1;
    public double w2;
    public double raised_w;
    public double downsampling;
    public String maskFilename;
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

        options.addOption(COMMAND_PARAMETERS.DIRECTORY, true, "");
        options.addOption(COMMAND_PARAMETERS.PORT, true, "");
        options.addOption(COMMAND_PARAMETERS.SINGLE_SELECTION, false, "");
        options.addOption(COMMAND_PARAMETERS.SELECTION_TYPE, true, "");
        options.addOption(COMMAND_PARAMETERS.FILTER, true, "");

        options.addOption(COMMAND_PARAMETERS.OPENING_MODE, true, "");
        options.addOption(COMMAND_PARAMETERS.POLL, false, "");
        options.addOption(COMMAND_PARAMETERS.ZOOM, true, "");
        options.addOption(COMMAND_PARAMETERS.RENDER_IMAGES, false, "");
        options.addOption(COMMAND_PARAMETERS.TABLE_ROWS, true, "");
        options.addOption(COMMAND_PARAMETERS.TABLE_COLUMNS, true, "");

        options.addOption(COMMAND_PARAMETERS.BAD_PIXELS_FACTOR, true, "");

        options.addOption(COMMAND_PARAMETERS.W1, true, "");
        options.addOption(COMMAND_PARAMETERS.W2, true, "");
        options.addOption(COMMAND_PARAMETERS.RAISED_W, true, "");

        options.addOption(COMMAND_PARAMETERS.DOWNSAMPLING, true, "");

        options.addOption(COMMAND_PARAMETERS.MASKFILENAME, true, "");

        options.addOption(COMMAND_PARAMETERS.INPUT_VECTORSFILE, true, "");
        options.addOption(COMMAND_PARAMETERS.INPUT_CLASSESFILE, true, "");
        options.addOption(COMMAND_PARAMETERS.INPUT_DATAFILE, true, "");


        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, args);

            if (cmdLine.hasOption(COMMAND_PARAMETERS.FILE)) {
                files = cmdLine.getOptionValues(COMMAND_PARAMETERS.FILE);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.DIRECTORY)) {
                directory = cmdLine.getOptionValue(COMMAND_PARAMETERS.DIRECTORY);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.PORT)) {
                port = Integer.parseInt(cmdLine.getOptionValue(COMMAND_PARAMETERS.PORT));
            }

            singleSelection = cmdLine.hasOption(COMMAND_PARAMETERS.SINGLE_SELECTION);

            if (cmdLine.hasOption(COMMAND_PARAMETERS.SELECTION_TYPE)) {
                selectionType = cmdLine.getOptionValue(COMMAND_PARAMETERS.SELECTION_TYPE);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.FILTER)) {
                String filters[] = cmdLine.getOptionValue(COMMAND_PARAMETERS.FILTER).split(COMMAND_PARAMETERS.FILTERS_SEPARATOR);

                filter = "";
                for (int i = 0; i < filters.length; i++) {
                    filter += filters[i];
                    if (i < filters.length - 1) {
                        filter += " ";
                    }
                }
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

            if (cmdLine.hasOption(COMMAND_PARAMETERS.BAD_PIXELS_FACTOR)) {
                bad_pixels_factor = Double.parseDouble(cmdLine.getOptionValue(COMMAND_PARAMETERS.BAD_PIXELS_FACTOR));
                System.out.println("bad_pixels.-.. "+bad_pixels_factor);
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.W1)) {
                w1 = Double.parseDouble(cmdLine.getOptionValue(COMMAND_PARAMETERS.W1));
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.W2)) {
                w2 = Double.parseDouble(cmdLine.getOptionValue(COMMAND_PARAMETERS.W2));
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.RAISED_W)) {
                raised_w = Double.parseDouble(cmdLine.getOptionValue(COMMAND_PARAMETERS.RAISED_W));
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.DOWNSAMPLING)) {
                downsampling = Double.parseDouble(cmdLine.getOptionValue(COMMAND_PARAMETERS.DOWNSAMPLING));
            }

            if (cmdLine.hasOption(COMMAND_PARAMETERS.MASKFILENAME)) {
                maskFilename = cmdLine.getOptionValue(COMMAND_PARAMETERS.MASKFILENAME);
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
        	//print error
        }
    }
}
