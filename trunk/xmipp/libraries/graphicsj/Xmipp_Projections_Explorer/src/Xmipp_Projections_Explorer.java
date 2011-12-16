
import window.JFrameLoad;
import constants.LABELS;
import explorer.ProjectionsExplorer;
import ij.IJ;
import ij.Macro;
import ij.plugin.PlugIn;
import java.io.File;
import java.io.IOException;
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
public class Xmipp_Projections_Explorer implements PlugIn {

    final static String COMMAND_OPTION_INPUT = "i";
    final static String COMMAND_OPTION_EULER_ANGLES = "angles";
    final static int INDEX_VOLUME = 0;
    final static int INDEX_EULER_ANGLES = 1;
    boolean use_sphere;
//    Sphere sphere;
//    Image3DUniverse universeVolume, universeSphere;
//    final static int UNIVERSE_W = 400, UNIVERSE_H = 400;
//    ImagePlus volumeIP, sphereIP;
//    ImageGeneric xmippVolume;  // Volume for Xmipp library.
//    boolean dispatched = false;
//    static ProjectionWindow projectionWindow;
//    JFrameImagesTable frameImagesTable;

    public static void main(String args[]) {
        //new ImageJ();
        String fileVolume = "/home/jvega/geoPhantom_128.vol";//phantom_128.vol";
        String fileEulerAngles = null;//"/home/jvega/temp/angles.sel";

        ProjectionsExplorer pe = new ProjectionsExplorer();
        pe.run(fileVolume, fileEulerAngles);
    }

    public void run(String string) {
        String fileVolume = null;
        String fileEulerAngles = null;
        boolean runGui = false;

        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "string" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            String argsList[] = processArgs(Macro.getOptions().trim());

            if (argsList[INDEX_VOLUME] != null) {
                try {
                    fileVolume = (new File(argsList[INDEX_VOLUME])).getAbsoluteFile().getCanonicalPath();
                } catch (IOException ioex) {
                    System.out.println("Error loading " + argsList[INDEX_VOLUME]);
                }
            } else {
                fileVolume = null;
            }

            if (argsList[INDEX_EULER_ANGLES] != null) {
                try {
                    fileEulerAngles = (new File(argsList[INDEX_EULER_ANGLES])).getAbsoluteFile().getCanonicalPath();
                } catch (IOException ioex) {
                    System.out.println("Error loading " + argsList[INDEX_EULER_ANGLES]);
                }
            } else {
                fileEulerAngles = null;
            }

            if (fileVolume == null) {
                runGui = true;
            } else {
                runGui = false;
            }
        } else {    // From menu.
            runGui = true;
        }

        if (runGui) {
            JFrameLoad frameLoad = new JFrameLoad();

            frameLoad.setLocationRelativeTo(null);
            frameLoad.setVisible(true);
            if (!frameLoad.cancelled) {
                fileVolume = frameLoad.getVolumeFile();
                fileEulerAngles = frameLoad.getEulerAnglesFile();
            }
        }

        // Volume: required, Angles: optional.
        if (fileVolume != null) {
            if (fileEulerAngles != null) {
                use_sphere = true;
                File f = new File(fileEulerAngles);
                if (!f.exists()) {
                    use_sphere = false;
                    IJ.error(LABELS.MESSAGE_ERROR_FNF(fileVolume));
                }
            }

            ProjectionsExplorer pe = new ProjectionsExplorer();
            pe.run(fileVolume, fileEulerAngles);
        } else {
            IJ.error(LABELS.MESSAGE_ERROR_FNF(fileVolume));
        }
    }

    static String[] processArgs(String args) {
        String argsList[] = args.split(" ");

        String parameters[] = {null, null};

        Options options = new Options();
        options.addOption(COMMAND_OPTION_INPUT, true, "Volume file");
        options.addOption(COMMAND_OPTION_EULER_ANGLES, true, "Euler angles file");

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            // Volume file.
            if (cmdLine.hasOption(COMMAND_OPTION_INPUT)) {
                parameters[INDEX_VOLUME] = cmdLine.getOptionValue(COMMAND_OPTION_INPUT);
            }

            // Angles file.
            if (cmdLine.hasOption(COMMAND_OPTION_EULER_ANGLES)) {
                parameters[INDEX_EULER_ANGLES] = cmdLine.getOptionValue(COMMAND_OPTION_EULER_ANGLES);
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }

        return parameters;
    }
}
