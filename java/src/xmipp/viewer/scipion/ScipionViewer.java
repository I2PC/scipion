/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.scipion;

import javax.swing.SwingUtilities;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.Param;
import xmipp.utils.XmippDialog;
import xmipp.viewer.Viewer;
import xmipp.viewer.windows.GalleryJFrame;
import static xmipp.viewer.windows.ImagesWindowFactory.openFileAsDefault;
import static xmipp.viewer.windows.ImagesWindowFactory.openFileAsImage;
import static xmipp.viewer.windows.ImagesWindowFactory.openMetadata;

/**
 *
 * @author airen
 */
public class ScipionViewer extends Viewer {

    public static void main(String args[]) {
        // Schedule a job for the event dispatch thread:
        // creating and showing this application's GUI.
        final String[] myargs = args;

        Param parameters = new Param(myargs);
        try {
            if (parameters.debug) {
                DEBUG.enableDebug(true);
            }

            if (parameters.files != null) {
                // DEBUG.enableDebug(true);
                final String[] files = parameters.files;

                for (int i = 0; i < files.length; i++) 
                    openFile(files[i], parameters);
                
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public static void openFile(String filename, Param parameters) {
        try {
            if (Filename.isMetadata(filename)) {
                if (parameters.mode.equalsIgnoreCase(Param.OPENING_MODE_IMAGE)) {
                    openFileAsImage(null, filename, parameters);
                } else {
                    parameters.mode = Param.OPENING_MODE_GALLERY;
                    openScipionGalleryJFrame(filename, parameters);
                }
            } else {
                ImageGeneric img = new ImageGeneric(filename);

                if (img.isSingleImage()) {
                    openFileAsImage(null, filename, parameters);
                } else if (img.isStackOrVolume()) {
                    if (parameters.mode
                            .equalsIgnoreCase(Param.OPENING_MODE_IMAGE)) {
                        openFileAsImage(null, filename, parameters);
                    } else {
                        parameters.mode = Param.OPENING_MODE_GALLERY;
                        openScipionGalleryJFrame(filename, parameters);
                    }
                }
            }
        } catch (Exception e) {
            XmippDialog.showError(null, String.format(
                    "Couldn't open file: '%s'\nError: %s", filename,
                    e.getMessage()));
            DEBUG.printException(e);
        }
    }

    public static void openScipionGalleryJFrame(final String filename, final Param parameters) throws Exception {
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {

                new ScipionGalleryJFrame(filename, new MetaData(filename), parameters);
            }
        });

    }

}
