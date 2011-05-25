/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package browser.windows;

import browser.imageitems.ImageConverter;
import ij.IJ;
import ij.ImagePlus;
import browser.table.JFrameImagesTable;
import browser.table.micrographs.ctf.CTFImageWindow;
import browser.table.micrographs.ctf.tasks.TasksEngine;
import ij.gui.ImageWindow;
import ij.gui.Toolbar;
import java.io.File;
import xmipp.Filename;
import xmipp.ImageDouble;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesWindowFactory {

    private final static String TEMPDIR_PATH = System.getProperty("java.io.tmpdir");

    public static void openFilesDefault(String files[], boolean poll) {
        for (int i = 0; i < files.length; i++) {
            openDefault(files[i], poll);
        }
    }

    public static void openDefault(String filename) {
        openDefault(filename, false);
    }

    public static void openDefault(String filename, boolean poll) {
        if (Filename.isSingleImage(filename)) {
            openFileAsImage(filename, poll);
        } else if (Filename.isStackOrVolume(filename)) {
            openFileAsTable(filename);
        } else if (Filename.isMetadata(filename)) {
            openFileAsTable(filename);
        } else {
            openFileAsImage(filename, poll);
        }
    }

    public static void openFilesAsImages(String filenames[], boolean poll) {
        for (int i = 0; i < filenames.length; i++) {
            String filename = Filename.getFilename(filenames[i]);
            long nimage = Filename.getNimage(filenames[i]);

            System.out.println(" *** Opening: " + filename + " / nimage: " + nimage);

            openFileAsImage(filename, poll);
        }
    }

    public static void openFileAsImage(String path) {
        openFileAsImage(path, false);
    }

    public static void openFileAsImage(String path, boolean poll) {
        File f = new File(path);

        if (f.exists()) {
            ImagePlus ip = null;

            if (Filename.isMetadata(path)) {
                ip = openMetaDataAsImage_(path);
            } else if (Filename.isXmippType(path)) {
                ip = openFileAsImage_(path);
            } else {
                ip = IJ.openImage(path);
            }

            openXmippImageWindow(ip, poll);
        } else {
            IJ.error("File not found: " + path);
        }
    }

    private static ImagePlus openFileAsImage_(String path) {
        ImagePlus ip = null;

        try {
            ImageDouble image = new ImageDouble(path);
            ip = ImageConverter.convertToImagej(image, path);
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        return ip;
    }

    private static ImagePlus openMetaDataAsImage_(String path) {
        ImagePlus ip = null;

        try {
            MetaData md = new MetaData(path);
            ip = ImageConverter.convertToImagej(md);
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        return ip;
    }

    private static ImageWindow openXmippImageWindow(ImagePlus imp, boolean poll) {
        ImageWindow iw = null;

        if (imp != null) {
            if (imp.getStackSize() > 1) {
                iw = new StackWindowOperations(imp, poll);
            } else {
                iw = new ImageWindowOperations(imp, poll);
            }
        }

        return iw;
    }

    public static void openFileAsTable(String filename) {
        File f = new File(filename);

        if (f.exists()) {
            JFrameImagesTable table = new JFrameImagesTable(filename);
            table.setVisible(true);
        } else {
            IJ.error("File not found: " + filename);
        }
    }

    public static void openFilesAsTable(String filenames[]) {
        openFilesAsTable(filenames, false);
    }

    public static void openFilesAsTable(String filenames[], boolean useSameTable) {
        if (useSameTable) {
            JFrameImagesTable table = new JFrameImagesTable(filenames);
            table.setVisible(true);
        } else {
            for (int i = 0; i < filenames.length; i++) {
                //JFrameImagesTable table = new JFrameImagesTable(filenames[i]);
                //table.setVisible(true);
                openFileAsTable(filenames[i]);
            }
        }
    }

    // Used by micrographs table, to load items marked as selected/unselected.
    public static void openTable(String filenames[], boolean enabled[]) {
        JFrameImagesTable table = new JFrameImagesTable(filenames, enabled);
        table.setVisible(true);
    }

    public static void captureFrame(ImagePlus ip) {
        openXmippImageWindow(ip, false);
    }

    public static ImageWindow openCTFImage(ImagePlus ip, String CTFfilename,
            String PSDfilename, TasksEngine tasksEngine,
            String MicrographFilename, int row) {
        IJ.setTool(Toolbar.FREEROI);

        return new CTFImageWindow(ip, CTFfilename, PSDfilename,
                tasksEngine, row);
    }
}
