

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import browser.imageitems.ImageConverter;
import browser.table.micrographs.JFrameMicrographs;
import ij.IJ;
import ij.Macro;
import ij.plugin.PlugIn;
import java.io.File;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import xmipp.ImageDouble;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class MicrographsBrowser implements PlugIn {

    private final static String COMMAND_OPTION_FILE = "i";

    public static void main_(String args[]) {
        try {
            String file = "/home/juanjo/Escritorio/images/image.spi";

            ImageDouble image = new ImageDouble();
            image.readPreview(file, 80, 80);
            ImageConverter.convertToImagej(image, file).show();
        } catch (Exception ex) {
            ex.printStackTrace();
            System.out.println("ERROR: " + ex.getMessage());
        }
    }

    public static void main(String args[]) {
        String filename = "/media/PENDRIVE/Ad5GLflagIIIa/all_micrographs_en.sel";

        openFile(filename);
    }

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

    private static void openFile(String fileName) {
        try {
            String path = "";
            if (!fileName.startsWith(File.separator)) {
                path = System.getProperty("user.dir") + File.separator;
            }

            path += fileName;

            MetaData md = new MetaData();
            md.read(path);

            JFrameMicrographs frameMicrographs = new JFrameMicrographs(md);
            frameMicrographs.setLocationRelativeTo(null);
            frameMicrographs.setVisible(true);
        } catch (Exception e) {
            e.printStackTrace();
            IJ.error(e.getMessage());
        }
    }

    private static String[] processArgs(String args) {
        String argsList[] = args.split(" ");
        Options options = new Options();

        options.addOption(COMMAND_OPTION_FILE, true, "file(s)");

        // It should be able to handle multiple files.
        options.getOption(COMMAND_OPTION_FILE).setOptionalArg(true);
        options.getOption(COMMAND_OPTION_FILE).setArgs(Integer.MAX_VALUE);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            if (cmdLine.hasOption(COMMAND_OPTION_FILE)) {
                return cmdLine.getOptionValues(COMMAND_OPTION_FILE);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return null;
    }
}
