/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.filebrowsers.model;

import ij.IJ;
import ij.io.Opener;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import javax.imageio.ImageIO;
import javax.imageio.stream.ImageInputStream;

/**
 *
 * @author Juanjo Vega
 */
public class FileBrowser {

    protected final static String BYTES = "bytes";
    protected final static String KBYTES = "KBytes";
    protected final static String MBYTES = "MBytes";
    protected final static String GBYTES = "GBytes";
    protected File currentDirectory;
    protected List<File> files;
    public final static String GIF_EXTS = ".gif";
    public final static String JPG_EXTS = ".jpg";
    public final static String JPEG_EXTS = ".jpeg";
    public final static String PNG_EXTS = ".png";
    public final static String BMP_EXTS = ".bmp";
//    public final static String IMAGE_EXTS[] = {GIF_EXTS, JPG_EXTS, JPEG_EXTS, PNG_EXTS, BMP_EXTS};

    public FileBrowser(String directory) {
        files = new ArrayList<File>();

        changeDirectory(new File(directory));
    }

    public String getCurrentDirectory() {
        try {
            return currentDirectory.getCanonicalPath();
        } catch (IOException ioex) {
            return null;
        }
    }

    public void refresh() {
        changeDirectory(currentDirectory);
    }

    public void goParent() {
        if (currentDirectory.getParentFile() != null) {
            changeDirectory(currentDirectory.getParentFile());
        }
    }

    public final void changeDirectory(File newDirectory) {
        try {
            currentDirectory = newDirectory.getCanonicalFile();

            File list[] = currentDirectory.listFiles();

            if (list == null) { // Avoiding exceptions when no files.
                files = Arrays.asList(new File[0]);
            } else {
                files = Arrays.asList(list);

                if (!files.isEmpty()) {
                    // Sorts files so folders will be at the top of list.
                    Collections.sort(files, new Comparator() {

                        public int compare(final Object o1, final Object o2) {
                            File a = (File) o1;
                            File b = (File) o2;

                            if (a.isDirectory() && !b.isDirectory()) {  // Dir vs. file
                                return -1;
                            } else if (!a.isDirectory() && b.isDirectory()) {   // file vs. dir
                                return 1;
                            } else {    // file vs. file
                                return a.compareTo(b);
                            }
                        }
                    });
                }
            }
        } catch (Exception ex) {
            IJ.error("Error", "Error changing to " + newDirectory.getAbsolutePath());
        }
    }

    public List<File> getFiles() {
        return files;
    }

    public static String getFileSizeString(double size) {
        String unit = BYTES;

        if (size > 1024) {  // bytes -> KB
            size /= 1024;
            unit = KBYTES;
            if (size > 1024) {  // KB -> MB
                size /= 1024;
                unit = MBYTES;
                if (size > 1024) {  // MB -> GB
                    size /= 1024;
                    unit = GBYTES;
                }
            }
        }

        DecimalFormat dc = new DecimalFormat(hasDecimalPart(size) ? ".##" : "");

        return dc.format(size) + " " + unit;
    }

    protected static boolean hasDecimalPart(double number) {
        return (number - (int) number) != 0;
    }

    public static boolean hasEnoughMemory(File file) {
        return file.length() < IJ.maxMemory();
    }

    public static boolean isFileSupported(File file) {
        try {
            final ImageInputStream iis = ImageIO.createImageInputStream(file);
            final Iterator readers = ImageIO.getImageReaders(iis);

            if (readers.hasNext()) {
                return true;
            } else {    // No image so...
                int type = (new Opener()).getFileType(file.getCanonicalPath());
                if (type != Opener.UNKNOWN) {   // ...let's see if it's a supported file type...
                    return true;
                }
            }
        } catch (Exception ex) {
        }

        return false;
    }

    protected static boolean isBasicImageType(File file) {
        if (isFileType(file, PNG_EXTS)
                || isFileType(file, GIF_EXTS)
                || isFileType(file, JPG_EXTS)
                || isFileType(file, JPEG_EXTS)
                || isFileType(file, BMP_EXTS)) {
            return true;
        }

        return false;
    }

    /**
     * Returns if given file is filetypes (by extension).
     * @param f
     * @param ext
     * @return
     */
    public static boolean isFileType(File f, String ext) {
        if (f.isDirectory()) {
            return false;
        }

        if (f.getName().toLowerCase().endsWith(ext.toLowerCase())) {
            return true;
        }

        return false;
    }
}
