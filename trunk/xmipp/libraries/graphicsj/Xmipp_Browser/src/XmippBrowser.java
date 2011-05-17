
import browser.JFrameBrowser;
import browser.LABELS;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import ij.Macro;
import ij.plugin.PlugIn;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import xmipp.MetaData;

/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippBrowser implements PlugIn {

    private final static String COMMAND_OPTION_DIR = "dir";
    private final static String COMMAND_OPTION_IMG = "img";
    private final static String COMMAND_OPTION_VOL = "vol";
    private final static String COMMAND_OPTION_SEL = "sel";
    private final static String COMMAND_OPTION_POLL = "poll";
    private String DIR;
    private String IMGS[];
    private String VOLS[];
    private String SELS[];
    private boolean POLL = false;
    protected JFrameBrowser frameBrowser;

    public void run(String string) {
        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "string" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            processArgs(Macro.getOptions().trim());
        } else {    // From menu.
            DIR = System.getProperty("user.dir");
        }

        if (IMGS != null) {
            ImagesWindowFactory.openImageFiles(IMGS, POLL);
        }

        if (VOLS != null) {
            for (int i = 0; i < VOLS.length; i++) {
                ImagesWindowFactory.openTable(VOLS[i]);//, h, w);
            }
        }

        if (SELS != null) {
            for (int i = 0; i < SELS.length; i++) {
                ImagesWindowFactory.openTable(SELS[i]);//, h, w);
            }
        }

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
        options.addOption(COMMAND_OPTION_VOL, true, "vol(s)");
        options.addOption(COMMAND_OPTION_SEL, true, "sel(s)");
        options.addOption(COMMAND_OPTION_POLL, false, "poll");

        // It should be able to handle multiple files.
        options.getOption(COMMAND_OPTION_IMG).setOptionalArg(true);
        options.getOption(COMMAND_OPTION_IMG).setArgs(Integer.MAX_VALUE);

        options.getOption(COMMAND_OPTION_VOL).setOptionalArg(true);
        options.getOption(COMMAND_OPTION_VOL).setArgs(Integer.MAX_VALUE);

        options.getOption(COMMAND_OPTION_SEL).setOptionalArg(true);
        options.getOption(COMMAND_OPTION_SEL).setArgs(Integer.MAX_VALUE);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            // Dir.
            if (cmdLine.hasOption(COMMAND_OPTION_DIR)) {
                DIR = cmdLine.getOptionValue(COMMAND_OPTION_DIR);
            }

            // Img(s).
            if (cmdLine.hasOption(COMMAND_OPTION_IMG)) {
                IMGS = cmdLine.getOptionValues(COMMAND_OPTION_IMG);

                String workdir = System.getProperty("user.dir");
                for (int i = 0; i < IMGS.length; i++) {
                    IMGS[i] = MetaData.fixPath(workdir, IMGS[i]);
                }
            }

            // Img(s).
            if (cmdLine.hasOption(COMMAND_OPTION_VOL)) {
                VOLS = cmdLine.getOptionValues(COMMAND_OPTION_VOL);
            }

            // Sel.
            if (cmdLine.hasOption(COMMAND_OPTION_SEL)) {
                SELS = cmdLine.getOptionValues(COMMAND_OPTION_SEL);
            }

            // Poll.
            if (cmdLine.hasOption(COMMAND_OPTION_POLL)) {
                POLL = true;
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }
}
