
import browser.JFrameBrowser;
import browser.LABELS;
import browser.table.JFrameImagesTable;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import ij.ImageJ;
import ij.Macro;
import ij.plugin.PlugIn;
import java.io.File;
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
    private final static String COMMAND_OPTION_IMG = "img";
    private final static String COMMAND_OPTION_SEL = "sel";
    private final static String COMMAND_OPTION_POLL = "poll";
    private String DIR;
    private String IMGS[];
    private String SEL[];
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
            DIR = System.getProperty("user.dir");
        }

        if (IMGS != null) {
            for (int i = 0; i < IMGS.length; i++) {
                if (!IMGS[i].isEmpty()) {
                    ImagesWindowFactory.openImageWindow(IMGS[i], POLL);
                }
            }
        }

        if (SEL != null) {
            for (int i = 0; i < SEL.length; i++) {
                if (!SEL[i].isEmpty()) {
//                    System.out.println(i + ": " + SEL[i]);

                    JFrameImagesTable selFileTable = new JFrameImagesTable();

                    File selfile = new File(SEL[i]);

                    String name = selfile.getName();
                    String parent = selfile.getParent();
                    if (parent == null) {
                        parent = System.getProperty("user.dir");
                    }

//                    System.out.println("Name: " + name);
//                    System.out.println("Parent: " + parent);

                    selFileTable.addSelFile(parent != null ? parent : DIR, name);

                    selFileTable.setVisible(true);
                }
            }
        }

//DIR = "/home/juanjo/Desktop/selfiles";
        if (DIR != null) {
            frameBrowser = new JFrameBrowser(LABELS.TITLE_MAIN_WINDOW, DIR);
            frameBrowser.setVisible(true);
        }
    }

    public void processArgs(String args) {
        String argsList[] = args.split(" ");
        Options options = new Options();

        options.addOption(COMMAND_OPTION_DIR, true, "directory");
        options.addOption(COMMAND_OPTION_IMG, true, "file(s)");
        options.addOption(COMMAND_OPTION_SEL, true, "sel(s)");
        options.addOption(COMMAND_OPTION_POLL, false, "poll");

        // It should be able to handle multiple files.
        options.getOption(COMMAND_OPTION_IMG).setOptionalArg(true);
        options.getOption(COMMAND_OPTION_IMG).setArgs(Integer.MAX_VALUE);

        options.getOption(COMMAND_OPTION_SEL).setOptionalArg(true);
        options.getOption(COMMAND_OPTION_SEL).setArgs(Integer.MAX_VALUE);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            // Dir.
            if (cmdLine.hasOption(COMMAND_OPTION_DIR)) {
                DIR = cmdLine.getOptionValue(COMMAND_OPTION_DIR);
            }

            // File(s).
            if (cmdLine.hasOption(COMMAND_OPTION_IMG)) {
                IMGS = cmdLine.getOptionValues(COMMAND_OPTION_IMG);
            }

            // Sel.
            if (cmdLine.hasOption(COMMAND_OPTION_SEL)) {
                SEL = cmdLine.getOptionValues(COMMAND_OPTION_SEL);
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
