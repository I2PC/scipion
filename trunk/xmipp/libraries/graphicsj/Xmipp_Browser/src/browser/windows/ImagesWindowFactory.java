/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package browser.windows;

import browser.DEBUG;
import browser.imageitems.ImageConverter;
import browser.imageitems.tableitems.AbstractTableImageItem;
import browser.table.models.AbstractXmippTableModel;
import browser.table.JFrameImagesTable;
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
import micrographs.ctf.tasks.TasksEngine;
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
        for (int i = 0; i < filenames.length; i++) {
            openFileAsDefault(filenames[i], poll);
        }
    }

    public static void openFileAsDefault(String filename) {
        openFileAsDefault(filename, false);
    }

    public static void openFileAsDefault(String filename, boolean poll) {
        if (Filename.isMetadata(filename)) {
            openFileAsTable(filename);
        } else {
            try {
                ImageDouble img = new ImageDouble();
                img.readHeader(filename);

                if (img.isSingleImage()) {
                    openFileAsImage(filename, poll);
                } else if (img.isStackOrVolume()) {
                    openFileAsTable(filename);
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

                imp = ImageConverter.convertToImagej(md);
            } else {
                ImageDouble id = new ImageDouble(path);
                imp = ImageConverter.convertToImagej(id, path);
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

    public static void openFileAsMetadata(String filename) {
        JFrameMetaData frameMetaData = new JFrameMetaData(filename);
        setConvenientSize(frameMetaData);
        frameMetaData.setLocationRelativeTo(null);
        frameMetaData.setVisible(true);
    }

    public static JFrameImagesTable openFileAsTable(String filename, int rows, int columns) {
        JFrameImagesTable table = new JFrameImagesTable(filename);
        setConvenientSize(table);

        table.setDimensions(rows, columns);
        table.setVisible(true);

        return table;
    }

    public static JFrameImagesTable openFileAsTable(String filename) {
        JFrameImagesTable table = new JFrameImagesTable(filename);
        setConvenientSize(table);

        table.setAutoAdjustColumns(true);
        table.setVisible(true);

        return table;
    }

    public static JFrameImagesTable openFilesAsTable(String filenames[]) {
        return openFilesAsTable(filenames, false);
    }

    public static JFrameImagesTable openFilesAsTable(String filenames[], boolean useSameTable) {
        JFrameImagesTable table = null;

        if (useSameTable) {
            table = new JFrameImagesTable(filenames);
            table.setAutoAdjustColumns(true);
            table.setVisible(true);

        } else {
            for (int i = 0; i < filenames.length; i++) {
                table = openFileAsTable(filenames[i]);
            }
        }

        return table;
    }

    // Used by micrographs table, to load items marked as selected/unselected.
    public static void openTable(String filenames[], boolean enabled[]) {
        JFrameImagesTable table = new JFrameImagesTable(filenames, enabled);
        table.setAutoAdjustColumns(true);
        table.setVisible(true);
    }

    public static void captureFrame(ImagePlus ip) {
        openXmippImageWindow(ip, false);
    }

    public static void openTableAs3D(AbstractXmippTableModel tableModel) {
        try {
            ArrayList<AbstractTableImageItem> items = tableModel.getAllItems();
            ImagePlus ip = ImageConverter.convertToImagePlus(items);
            ip.setTitle(tableModel.getFilename());

            openImagePlusAs3D(ip);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            ex.printStackTrace();
        }
    }

    public static void openTableAsImagePlus(AbstractXmippTableModel tableModel) {
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

                ArrayList<AbstractTableImageItem> items = tableModel.getAllItems();
                ImagePlus imp = ImageConverter.convertToImagePlus(items);
                IJ.run(imp, "Xmipp writer", "save=" + tempFile.getAbsolutePath());

//                System.err.println(" >>> TMP Saved at: " + file.getAbsolutePath());

                imp.setTitle(tempFile.getName());

                captureFrame(imp);
            }
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            ex.printStackTrace();
        }
    }

    public static void openImagePlusAsTable(ImagePlus imp) {
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

            openFileAsTable(file.getAbsolutePath());
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            ex.printStackTrace();
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

    public static void openFileAsText(String filename, Component parent) {
        JFrameTextFile frameCTF = new JFrameTextFile(filename);
        frameCTF.setLocationRelativeTo(parent);
        frameCTF.setVisible(true);
    }

    public static void openCTFView(ImagePlus imp, String CTFFilename, String PSDFilename) {
        CTFProfileWindow ctfView = new CTFProfileWindow(imp, CTFFilename, PSDFilename);
        ctfView.setVisible(true);
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
