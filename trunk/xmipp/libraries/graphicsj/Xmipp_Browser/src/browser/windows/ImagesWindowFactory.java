/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package browser.windows;

import browser.DEBUG;
import browser.imageitems.ImageConverter;
import browser.imageitems.tableitems.AbstractGalleryImageItem;
import browser.gallery.JFrameGallery;
import browser.gallery.models.AbstractXmippTableModel;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageWindow;
import ij.gui.Toolbar;
import ij.io.FileInfo;
import ij.process.StackConverter;
import ij3d.Content;
import ij3d.Image3DUniverse;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Toolkit;
import java.io.File;
import java.util.ArrayList;
import javax.swing.JFrame;
import metadata.JFrameMetaData;
import micrographs.JFrameMicrographs;
import micrographs.ctf.CTFRecalculateImageWindow;
import micrographs.CTFProfileWindow;
import micrographs.FSCWindow;
import micrographs.ctf.tasks.TasksEngine;
import rotspectra.JFrameRotSpectra;
import xmipp.Filename;
import xmipp.ImageDouble;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesWindowFactory {

    private final static int UNIVERSE_W = 400, UNIVERSE_H = 400;
//    private final static String TEMPDIR_PATH = System.getProperty("java.io.tmpdir");

    public static void openFilesAsDefault(String filenames[], boolean poll) {
        openFilesAsDefault(filenames, poll, -1, -1);
    }

    public static void openFilesAsDefault(String filenames[], int rows, int columns) {
        openFilesAsDefault(filenames, false, rows, columns);
    }

    public static void openFilesAsDefault(String filenames[], boolean poll, int rows, int columns) {
        for (int i = 0; i < filenames.length; i++) {
            openFileAsDefault(filenames[i], poll, rows, columns);
        }
    }

    public static void openFileAsDefault(String filenames) {
        openFileAsDefault(filenames, -1, -1);
    }

    public static void openFileAsDefault(String filename, int rows, int columns) {
        openFileAsDefault(filename, false, rows, columns);
    }

    public static void openFileAsDefault(String filename, boolean poll, int rows, int columns) {
        if (Filename.isMetadata(filename)) {
            openFileAsGallery(filename, rows, columns);
        } else {
            try {
                ImageDouble img = new ImageDouble();
                img.readHeader(filename);

                if (img.isSingleImage()) {
                    openFileAsImage(filename, poll);
                } else if (img.isStackOrVolume()) {
                    openFileAsGallery(filename, rows, columns);
                } else {
                    openFileAsImage(filename, poll);
                }
            } catch (Exception e) {
                IJ.open(filename);
                //IJ.error(e.getMessage());
            }
        }
    }

    public static void openFilesAsImages(String filenames[], boolean poll) {
        for (int i = 0; i < filenames.length; i++) {
            String filename = Filename.getFilename(filenames[i]);
            long nimage = Filename.getNimage(filenames[i]);

            DEBUG.printMessage(" *** Opening: " + filename + " / nimage: " + nimage);

            openFileAsImage(filenames[i], poll);
        }
    }

    public static void openFileAsImage(String path) {
        openFileAsImage(path, false);
    }

    public static void openFileAsImage(String path, boolean poll) {
        try {
            ImagePlus imp;

            if (Filename.isMetadata(path)) {
                MetaData md = new MetaData(path);

                imp = ImageConverter.convertToImageJ(md);
            } else {
                ImageDouble id = new ImageDouble(path);
                imp = ImageConverter.convertToImageJ(id, path);
            }

            // Normalize image stack.
            ImageConverter.normalizeStack(imp);

            openXmippImageWindow(imp, poll);
        } catch (Exception ex) {
            IJ.error(ex.getMessage() + ": " + path);
            DEBUG.printException(ex);
        }
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

    public static void openFilesAsMetadata(String filenames[]) {
        openFilesAsMetadata(filenames, false);
    }

    public static void openFilesAsMetadata(String filenames[], boolean render) {
        for (int i = 0; i < filenames.length; i++) {
            openFileAsMetadata(filenames[i], render);
        }
    }

    public static void openFileAsMetadata(String filename) {
        openFileAsMetadata(filename, false);
    }

    public static void openFileAsMetadata(String filename, boolean render) {
        JFrameMetaData frameMetaData = new JFrameMetaData(filename);
        frameMetaData.setRenderImages(render);
        setConvenientSize(frameMetaData);
        frameMetaData.setLocationRelativeTo(null);
        frameMetaData.setVisible(true);
    }

    public static void openFileAsGallery(String filename) {
        openFileAsGallery(filename, -1, -1);
    }

    public static JFrameGallery openFileAsGallery(String filename, int rows, int columns) {
        JFrameGallery gallery = new JFrameGallery(filename);
        setConvenientSize(gallery);

        if (rows < 0 && columns < 0) {
            gallery.setAutoAdjustColumns(true);
        } else {
            gallery.setDimensions(rows, columns);
        }

        gallery.setVisible(true);

        return gallery;
    }

    public static JFrameGallery openFilesAsGallery(String filenames[], boolean useSameTable) {
        return openFilesAsGallery(filenames, useSameTable, -1, -1);
    }

    public static JFrameGallery openFilesAsGallery(String filenames[], int rows, int columns) {
        return openFilesAsGallery(filenames, false, rows, columns);
    }

    public static JFrameGallery openFilesAsGallery(String filenames[],
            boolean useSameTable, int rows, int columns) {
        JFrameGallery gallery = null;

        if (useSameTable) {
            gallery = new JFrameGallery(filenames);
            setConvenientSize(gallery);

            if (rows < 1 && columns < 1) {
                gallery.setAutoAdjustColumns(true);
            } else {
                gallery.setDimensions(rows, columns);
            }
            gallery.setVisible(true);
        } else {
            for (int i = 0; i < filenames.length; i++) {
                gallery = openFileAsGallery(filenames[i], rows, columns);
            }
        }

        return gallery;
    }

    // Used by micrographs table, to load items marked as selected/unselected.
    public static void openGallery(String filenames[], boolean enabled[]) {
        JFrameGallery gallery = new JFrameGallery(filenames, enabled);
        gallery.setAutoAdjustColumns(true);
        gallery.setVisible(true);
    }

    public static void captureFrame(ImagePlus ip) {
        openXmippImageWindow(ip, false);
    }

    public static void openGalleryAs3D(AbstractXmippTableModel tableModel) {
        try {
            ArrayList<AbstractGalleryImageItem> items = tableModel.getAllItems();
            ImagePlus ip = ImageConverter.convertToImageJ(items);
            ip.setTitle(tableModel.getFilename());

            openImagePlusAs3D(ip);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            DEBUG.printException(ex);
        }
    }

    public static void openGalleryAsImagePlus(AbstractXmippTableModel tableModel) {
        try {
            String path = tableModel.getFilename();

            // If there is an associated filename, uses it...
            File file = new File(path);
            if (file.exists()) {
//                System.err.println(" +++ EXISTS");
                openFileAsImage(path);
            } else {
//                System.err.println(" !!! EXISTS");
                // ...otherwise, stores it in a temporary file.
                File tempFile = File.createTempFile("tableToStack_", ".stk");
                tempFile.deleteOnExit();

                ArrayList<AbstractGalleryImageItem> items = tableModel.getAllItems();
                ImagePlus imp = ImageConverter.convertToImageJ(items);
                IJ.run(imp, "Xmipp writer", "save=" + tempFile.getAbsolutePath());

//                System.err.println(" >>> TMP Saved at: " + file.getAbsolutePath());

                imp.setTitle(tempFile.getName());

                captureFrame(imp);
            }
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            DEBUG.printException(ex);
        }
    }

    public static void openImagePlusAsGallery(ImagePlus imp) {
        try {
            FileInfo fi = imp.getOriginalFileInfo();
//            System.out.println(" +++ FileInfo: " + fi);

            // If path exists, uses it...
            File file = null;
            if (fi != null && !fi.fileName.trim().isEmpty() && !fi.directory.trim().isEmpty()) {
                file = new File(fi.directory + File.separator + fi.fileName);
            }

            if (file == null || !file.exists()) {   // ...otherwise, stores it in a temporary file.
//                System.err.println(" !!! EXISTS");
                file = File.createTempFile("stackToTable_", ".stk");
                file.deleteOnExit();
                IJ.run(imp, "Xmipp writer", "save=" + file.getAbsolutePath());

//                System.err.println(" >>> TMP Saved at: " + file.getAbsolutePath());
//            } else {
//                System.err.println(" +++ EXISTS");
            }

            openFileAsGallery(file.getAbsolutePath(), -1, -1);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            DEBUG.printException(ex);
        }
    }

    public static void openImagePlusAs3D(ImagePlus ip) {
        Image3DUniverse universe = new Image3DUniverse(UNIVERSE_W, UNIVERSE_H);

        // Adds the sphere image plus to universe.
        new StackConverter(ip).convertToRGB();
        Content c = universe.addVoltex(ip);
        c.displayAs(Content.VOLUME);

        universe.show();    // Shows...
    }

    public static ImageWindow openCTFImage(ImagePlus ip, String CTFfilename,
            String PSDfilename, TasksEngine tasksEngine,
            String MicrographFilename, int row) {
        IJ.setTool(Toolbar.FREEROI);

        return new CTFRecalculateImageWindow(ip, CTFfilename, PSDfilename,
                tasksEngine, row);
    }

    public static void openMicrograph(String filename) {
        File f = new File(filename);

        if (f.exists()) {
            JFrameMicrographs frame = new JFrameMicrographs(filename);
            frame.setVisible(true);
        } else {
            IJ.error("File is missing", filename + " not found.");
        }
    }

    public static void openRotSpectrasWindow(String filenameVectors,
            String filenameClasses, String filenameData) {
        JFrameRotSpectra frame = new JFrameRotSpectra(filenameVectors, filenameClasses, filenameData);
        ImagesWindowFactory.setConvenientSize(frame);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    public static void openFileAsText(String filename, Component parent) {
        JFrameTextFile frameCTF = new JFrameTextFile(filename);
        frameCTF.setLocationRelativeTo(parent);
        frameCTF.setVisible(true);
    }

    public static void openCTFWindow(ImagePlus imp, String CTFFilename, String PSDFilename) {
        CTFProfileWindow ctfView = new CTFProfileWindow(imp, CTFFilename, PSDFilename);
        ctfView.setVisible(true);
    }

    public static void openFSCWindow(String filename) {
        FSCWindow frame = new FSCWindow(filename);
        frame.setVisible(true);
    }

    public static String getSortTitle(String title, int width, FontMetrics fontMetrics) {
        String sort = title;
        int strlenght = fontMetrics.stringWidth(sort);
        int index = 0;

        while (strlenght > width) {
            index++;
            sort = "..." + title.substring(index);
            strlenght = fontMetrics.stringWidth(sort);
        }

        return sort;
    }

    public static void setConvenientSize(JFrame frame) {
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        int w = screenSize.width * 2 / 3;
        int h = screenSize.height * 2 / 3;

        frame.setSize(w, h);
        frame.setLocationRelativeTo(null);
    }
}
