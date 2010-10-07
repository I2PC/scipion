
import browser.JFrameBrowser;
import browser.LABELS;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import ij.ImageJ;
import ij.Macro;
import ij.WindowManager;
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
public class XmippBrowser implements PlugIn {

    private final static String COMMAND_OPTION_DIR = "dir";
    private final static String COMMAND_OPTION_FILE = "file";
    private final static String COMMAND_OPTION_POLL = "poll";
    private String DIR;
    private String FILES[];
    private boolean POLL = false;
    protected JFrameBrowser frameBrowser;

    public static void main(String args[]) {
        new ImageJ();

        (new XmippBrowser()).run("");
    }

    public XmippBrowser() {
        //   frameBrowser = new JFrameBrowser(LABELS.TITLE_MAIN_WINDOW, work_dir);
    }

    public void run(String string) {
        if (IJ.isMacro()) { // From macro.
            // "string" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            processArgs(Macro.getOptions().trim());
        } else {    // From menu.
        }

        String work_dir = null;
        if (DIR != null) {
            work_dir = DIR;
        } else if (FILES == null) {
            work_dir = System.getProperty("user.dir");
        }

        if (FILES != null) {
            for (int i = 0; i < FILES.length; i++) {
                if (!FILES[i].isEmpty()) {
                    System.out.println("Opening: " + FILES[i]);
                    ImagesWindowFactory.openImageWindow(FILES[i], POLL);
                }
            }
        }

        if (work_dir != null) {
            frameBrowser = new JFrameBrowser(LABELS.TITLE_MAIN_WINDOW, work_dir);
            frameBrowser.setVisible(true);
        }
    }

    public void processArgs(String args) {
        String argsList[] = args.split(" ");
        Options options = new Options();

        options.addOption(COMMAND_OPTION_DIR, true, "directory");
        options.addOption(COMMAND_OPTION_FILE, true, "file(s)");
        options.addOption(COMMAND_OPTION_POLL, false, "poll");

        // It should be able to handle multiple files.
        options.getOption(COMMAND_OPTION_FILE).setOptionalArg(true);
        options.getOption(COMMAND_OPTION_FILE).setArgs(Integer.MAX_VALUE);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            // Dir.
            if (cmdLine.hasOption(COMMAND_OPTION_DIR)) {
                DIR = cmdLine.getOptionValue(COMMAND_OPTION_DIR);
            }

            // File(s).
            if (cmdLine.hasOption(COMMAND_OPTION_FILE)) {
                FILES = cmdLine.getOptionValues(COMMAND_OPTION_FILE);
            }

            // Label.
            if (cmdLine.hasOption(COMMAND_OPTION_POLL)) {
                POLL = true;
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
