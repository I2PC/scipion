

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import browser.COMMAND_PARAMETERS;
import metadata.JFrameMetaData;
import ij.IJ;
import ij.Macro;
import ij.plugin.PlugIn;
import java.io.File;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

/**
 *
 * @author Juanjo Vega
 */
public class XmippMetaDataViewer implements PlugIn {

    String INPUT[];
    boolean RENDER_IMAGES;

    public void run(String args) {
        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "args" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            processArgs(Macro.getOptions().trim());

            if (INPUT != null) {
                for (int i = 0; i < INPUT.length; i++) {
                    if (!INPUT[i].isEmpty()) {
                        final String fileName = INPUT[i];
                        java.awt.EventQueue.invokeLater(new Runnable() {

                            public void run() {
                                openFile(fileName, RENDER_IMAGES);
                            }
                        });
                    }
                }
            }
        }
    }

    public static void openFile(String fileName, boolean render_images) {
        try {
            String path = "";
            if (!fileName.startsWith(File.separator)) {
                path = System.getProperty("user.dir") + File.separator;
            }

            path += fileName;

            JFrameMetaData frameMetaData = new JFrameMetaData(path);
            frameMetaData.setRenderImages(render_images);
            frameMetaData.setLocationRelativeTo(null);
            frameMetaData.setVisible(true);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    void processArgs(String args) {
        String argsList[] = args.split(" ");
        Options options = new Options();

        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE, true, "");
        options.addOption(COMMAND_PARAMETERS.OPTION_RENDER_IMAGES, false, "");

        // It should be able to handle multiple files.
        options.getOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE).setOptionalArg(true);
        options.getOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE).setArgs(Integer.MAX_VALUE);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            // Input.
            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE)) {
                INPUT = cmdLine.getOptionValues(COMMAND_PARAMETERS.OPTION_INPUT_FILE);
            }

            // Input.
            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_RENDER_IMAGES)) {
                RENDER_IMAGES = true;
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }
}
