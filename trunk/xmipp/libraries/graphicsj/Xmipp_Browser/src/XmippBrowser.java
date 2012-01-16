
import browser.DEBUG;
import ij.IJ;
import ij.Macro;
import browser.commandline.Parameters;
import browser.commandline.COMMAND_PARAMETERS;
import browser.windows.ImagesWindowFactory;
import ij.plugin.PlugIn;
import java.util.LinkedList;
import xmipp.jni.Filename;

/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippBrowser implements PlugIn {

    public XmippBrowser() {
        DEBUG.enableDebug(false);
    }

    @Override
    public void run(String string) {
        Parameters parameters;

        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "string" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            parameters = new Parameters(Macro.getOptions().split(" "));
//            processArgs(Macro.getOptions().trim());
        } else {    // From menu.
            parameters = new Parameters();
        }

        if (parameters.files != null) {
            openFiles(parameters);
        }
    }

    static void openFiles(Parameters parameters) {
        if (parameters.mode.trim().toLowerCase().compareTo(COMMAND_PARAMETERS.OPENING_MODE_DEFAULT) == 0) {
            ImagesWindowFactory.openFilesAsDefault(parameters.files, parameters);
        } else if (parameters.mode.trim().toLowerCase().compareTo(COMMAND_PARAMETERS.OPENING_MODE_IMAGE) == 0) {
            ImagesWindowFactory.openFilesAsImages(parameters.files, parameters);
        } else if (parameters.mode.trim().toLowerCase().compareTo(COMMAND_PARAMETERS.OPENING_MODE_GALLERY) == 0) {
            // Separates single images (index 0) from the rest of files (index 1).
            String filesList[][] = getSeparatedFiles(parameters.files);

            // Single images are opened in the same table...
            if (filesList[0] != null) {
                ImagesWindowFactory.openFilesAsGallery(filesList[0], true, parameters);
            }

            // ...while rest of files are opened in separated ones.
            if (filesList[1] != null) {
                ImagesWindowFactory.openFilesAsGallery(filesList[1], false, parameters);
            }
        } else if (parameters.mode.trim().toLowerCase().compareTo(COMMAND_PARAMETERS.OPENING_MODE_METADATA) == 0) {
            ImagesWindowFactory.openFilesAsMetadata(parameters.files, parameters);
        }
    }

    static String[][] getSeparatedFiles(String items[]) {
        LinkedList<String> images = new LinkedList<String>();
        LinkedList<String> files = new LinkedList<String>();
        LinkedList<String> storeTo;

        for (int i = 0; i < items.length; i++) {
            String filename = items[i];

            try {
                storeTo = Filename.isSingleImage(filename) ? images : files;
                storeTo.add(filename);
            } catch (Exception e) {
                IJ.error(e.getMessage());
            }
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
}
