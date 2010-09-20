
import browser.JFrameBrowser;
import browser.LABELS;
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
public class XmippBrowser implements PlugIn {

    private final static String COMMAND_OPTION_ZOOM = "zoom";
    private final static String COMMAND_OPTION_GO_TO_INDEX = "gotoindex";
    private final static String COMMAND_OPTION_GO_TO_LABEL = "gotolabel";
    private final static int INDEX_ZOOM = 0;
    private final static int INDEX_GO_TO_INDEX = 1;
    private final static int INDEX_GO_TO_LABEL = 2;
    protected JFrameBrowser frameBrowser;

    public static void main(String args[]) {
        new ImageJ();

        (new XmippBrowser()).run("");
    }

    public XmippBrowser() {
        //   frameBrowser = new JFrameBrowser(LABELS.TITLE_MAIN_WINDOW, work_dir);
    }

    public void run(String string) {
        String work_dir = System.getProperty("user.dir");
        /*        work_dir = "/home/juanjo/Desktop/images/12205" //"/home/juanjo/Desktop/SampleData/Imgs"
        //"/home/juanjo/Desktop/selfiles");
        //"C:\\Documents and Settings\\Juanjo\\Escritorio\\imgs");
        //System.getProperty("user.dir")
        ;*/

        if (IJ.isMacro()) { // From macro.
            // "string" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            //String argsList[] = processArgs(Macro.getOptions());

            work_dir = Macro.getOptions().trim();
            //System.out.println("[" + work_dir + "]");

//            fileVolume = argsList[INDEX_VOLUME];
//            fileEulerAngles = argsList[INDEX_EULER_ANGLES];
        } else {    // From menu.
        }

        frameBrowser = new JFrameBrowser(LABELS.TITLE_MAIN_WINDOW, work_dir);
        frameBrowser.setVisible(true);
    }

    public static String[] processArgs(String args) {
        String argsList[] = args.split(" ");

        Options options = new Options();
        options.addOption(COMMAND_OPTION_ZOOM, true, "Initial Zoom");
        options.addOption(COMMAND_OPTION_GO_TO_INDEX, true, "selected index");
        options.addOption(COMMAND_OPTION_GO_TO_LABEL, true, "selected label");

        String parameters[] = new String[options.getOptions().size()];

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            // Zoom.
            if (cmdLine.hasOption(COMMAND_OPTION_ZOOM)) {
                parameters[INDEX_ZOOM] = cmdLine.getOptionValue(COMMAND_OPTION_ZOOM);
            }

            // Index.
            if (cmdLine.hasOption(COMMAND_OPTION_GO_TO_INDEX)) {
                parameters[INDEX_GO_TO_INDEX] = cmdLine.getOptionValue(COMMAND_OPTION_GO_TO_INDEX);
            }

            // Label.
            if (cmdLine.hasOption(COMMAND_OPTION_GO_TO_LABEL)) {
                parameters[INDEX_GO_TO_LABEL] = cmdLine.getOptionValue(COMMAND_OPTION_GO_TO_LABEL);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return parameters;
    }
}
