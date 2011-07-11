
import browser.JFrameBrowser;
import browser.LABELS;
import browser.COMMAND_PARAMETERS;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import ij.Macro;
import ij.plugin.PlugIn;
import java.util.LinkedList;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import xmipp.Filename;

/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippBrowser implements PlugIn {

    private String DIR;
    private String INPUT[];
    private String MODE = COMMAND_PARAMETERS.MODE_DEFAULT;
    private boolean POLL = false;
    private final static int MODE_TYPE_DEFAULT = 0;
    private final static int MODE_TYPE_IMAGE = 1;
    private final static int MODE_TYPE_TABLE = 2;
    protected JFrameBrowser frameBrowser;

    public void run(String string) {
        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "string" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            processArgs(Macro.getOptions().trim());
        } else {    // From menu.
            DIR = System.getProperty("user.dir");
        }

        // @TODO Remove this.
        // -------------------------
//        DIR = null;
//        processArgs(string);
        // -------------------------

        if (INPUT != null) {
            openFiles(INPUT, MODE, POLL);
        }

        if (DIR != null) {
            frameBrowser = new JFrameBrowser(LABELS.TITLE_MAIN_WINDOW, DIR);
            frameBrowser.setVisible(true);
        }
    }

    private static void openFiles(String input[], String mode, boolean poll) {
        int m = 0;
        if (mode.trim().toLowerCase().compareTo(COMMAND_PARAMETERS.MODE_DEFAULT) == 0) {
            m = 0;
        } else if (mode.trim().toLowerCase().compareTo(COMMAND_PARAMETERS.MODE_IMAGE) == 0) {
            m = 1;
        } else if (mode.trim().toLowerCase().compareTo(COMMAND_PARAMETERS.MODE_TABLE) == 0) {
            m = 2;
        }

        if (input != null) {
            switch (m) {
                case MODE_TYPE_DEFAULT:
                    ImagesWindowFactory.openFilesAsDefault(input, poll);
                    break;
                case MODE_TYPE_IMAGE:
                    ImagesWindowFactory.openFilesAsImages(input, poll);
                    break;
                case MODE_TYPE_TABLE:
                    // Separates single images (index 0) from the rest of files (index 1).
                    String filesList[][] = getSeparatedFiles(input);

                    // Single images are opened in the same table...
                    if (filesList[0] != null) {
                        ImagesWindowFactory.openFilesAsTable(filesList[0], true);
                    }

                    // ...while rest of files are opened in separated ones.
                    if (filesList[1] != null) {
                        ImagesWindowFactory.openFilesAsTable(filesList[1]);
                    }
                    break;
            }
        }
    }

    private static String[][] getSeparatedFiles(String items[]) {
        LinkedList<String> images = new LinkedList<String>();
        LinkedList<String> files = new LinkedList<String>();
        LinkedList<String> storeTo;

        for (int i = 0; i < items.length; i++) {
            String filename = items[i];

            storeTo = Filename.isSingleImage(filename) ? images : files;
            storeTo.add(filename);
        }

        String result[][] = new String[2][];

        if (!images.isEmpty()) {
            result[0] = images.toArray(new String[images.size()]);
        }

        if (!files.isEmpty()) {
            result[1] = files.toArray(new String[files.size()]);
        }

        return result;
    }

    private void processArgs(String args) {
        String argsList[] = args.split(" ");
        Options options = new Options();

        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_DIR, true, COMMAND_PARAMETERS.OPTION_INPUT_DIR_DESCRIPTION);
        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE, true, COMMAND_PARAMETERS.OPTION_INPUT_FILE_DESCRIPTION);
        options.addOption(COMMAND_PARAMETERS.OPTION_MODE, true, COMMAND_PARAMETERS.OPTION_MODE_DESCRIPTION);
        options.addOption(COMMAND_PARAMETERS.OPTION_POLL, false, COMMAND_PARAMETERS.OPTION_POLL_DESCRIPTION);

        // It should be able to handle multiple files.
        options.getOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE).setOptionalArg(true);
        options.getOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE).setArgs(Integer.MAX_VALUE);

        options.getOption(COMMAND_PARAMETERS.OPTION_MODE).setOptionalArg(true);
        options.getOption(COMMAND_PARAMETERS.OPTION_MODE).setArgs(Integer.MAX_VALUE);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            // Dir.
            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_DIR)) {
                DIR = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_INPUT_DIR);
            }

            // Input.
            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE)) {
                String filenames[] = cmdLine.getOptionValues(COMMAND_PARAMETERS.OPTION_INPUT_FILE);
                LinkedList<String> list = new LinkedList<String>();

                String workdir = System.getProperty("user.dir");
                for (int i = 0; i < filenames.length; i++) {
                    if (!filenames[i].trim().isEmpty()) {
                        list.add(Filename.fixPath(workdir, filenames[i]));
                    }
                }

                INPUT = list.toArray(new String[list.size()]);
            }

            // Mode.
            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_MODE)) {
                MODE = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_MODE);
            }

            // Poll.
            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_POLL)) {
                POLL = true;
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }
}
