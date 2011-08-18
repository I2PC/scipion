

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import browser.COMMAND_PARAMETERS;
import ij.IJ;
import ij.Macro;
import ij.plugin.PlugIn;
import java.io.File;
import micrographs.JFrameMicrographs;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

/**
 *
 * @author Juanjo Vega
 */
public class XmippMicrographViewer implements PlugIn {

    public void run(String args) {
        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "args" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            String FILES[] = processArgs(Macro.getOptions().trim());

            if (FILES != null) {
                for (int i = 0; i < FILES.length; i++) {
                    if (!FILES[i].isEmpty()) {
                        final String fileName = FILES[i];
                        java.awt.EventQueue.invokeLater(new Runnable() {

                            public void run() {
                                openFile(fileName);
                            }
                        });
                    }
                }
            }
        }
    }

    public static void openFile(String fileName) {
        try {
            String path = "";
            if (!fileName.startsWith(File.separator)) {
                path = System.getProperty("user.dir") + File.separator;
            }

            path += fileName;

            JFrameMicrographs frameMicrographs = new JFrameMicrographs(path);
            frameMicrographs.setLocationRelativeTo(null);
            frameMicrographs.setVisible(true);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    private static String[] processArgs(String args) {
        String argsList[] = args.split(" ");
        Options options = new Options();

        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE, true, COMMAND_PARAMETERS.OPTION_INPUT_FILE_DESCRIPTION);

        // It should be able to handle multiple files.
        options.getOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE).setOptionalArg(true);
        options.getOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE).setArgs(Integer.MAX_VALUE);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE)) {
                return cmdLine.getOptionValues(COMMAND_PARAMETERS.OPTION_INPUT_FILE);
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }

        return null;
    }
}
