package xmipp;

import ij.plugin.*;
import ij.*;
import ij.io.OpenDialog;
import xmipp.ij.io.*;
import java.io.*;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: MESSAOUDI cedric
 * Date: 9 dec. 2008
 * Time: 11:06:03
 * To change this template use File | Settings | File Templates.
 */
/**
 * Modified Wed May 26, 2010
 * @author Juanjo Vega
 * Added support for new sel files format.
 */
public class Sel_Reader extends ImagePlus implements PlugIn {

    private final static String IMAGE = "image";
    private final static String ENABLED = "enabled";
    private final static int INDEX_IMAGE = 0;
    private final static int INDEX_ENABLED = 1;

    /**
     *  Main processing method for the SpiderReader_ object
     *
     *@param  arg  Description of the Parameter
     */
    public void run(String arg) {
        String directory = "";
        String fileName = arg;
        if ((arg == null) || (arg.compareTo("") == 0)) {
            // Choose a file since none specified
            OpenDialog od = new OpenDialog("Load sel File...", arg);
            fileName = od.getFileName();
            if (fileName == null) {
                return;
            }
            directory = od.getDirectory();
        } else {
            // we were sent a filename to open
            File dest = new File(arg);
            directory = dest.getParent();
            fileName = dest.getName();
        }

        // Load in the image
        ImagePlus imp = load(directory, fileName);
        if (imp == null) {
            return;
        }

        // Attach the Image Processor
        if (imp.getNSlices() > 1) {
            setStack(fileName, imp.getStack());
        } else {
            setProcessor(fileName, imp.getProcessor());
        }
        // Copy the scale info over
        copyScale(imp);

        // Show the image if it was selected by the file
        // chooser, don't if an argument was passed ie
        // some other ImageJ process called the plugin
        if (arg.equals("")) {
            show();
        }
    }

    /**
     *  Description of the Method
     *
     *@param  directory  Description of the Parameter
     *@param  fileName   Description of the Parameter
     *@return            Description of the Return Value
     */
    public static ImagePlus load(String directory, String fileName) {
        if (fileName == null || fileName.isEmpty()) {
            return null;
        }

        if (!directory.endsWith(File.separator)) {
            directory += File.separator;
        }
        ImageStack is = null;
        IJ.showStatus("Loading sel File: " + directory + fileName);
        try {
            RandomAccessFile f = new RandomAccessFile(directory + fileName, "r");
            Opener tr = new Opener();
            String line;
            boolean oldFormat;
            int indexes[] = null;

            // Skips header
            line = f.readLine();//; XMIPP_3 * column_format *
            if (line.toUpperCase().contains("XMIPP_3")) {
                oldFormat = false;

                f.readLine();// 2nd header line ";"

                // 3rd line contains parameters order.
                String metadata[] = f.readLine().split(" ");
                indexes = new int[metadata.length - 1];
                for (int i = 0; i < indexes.length; i++) {
                    indexes[i] = -1;
                }

                for (int i = 1; i < metadata.length; i++) { // 0 position is ";"
                    int index = getFieldIndex(metadata[i]);
                    if (index >= 0) {
                        indexes[index] = i - 1;
                    }
                }
            } else {  // Old format.
                oldFormat = true;
                f.seek(0);
            }

            while ((line = f.readLine()) != null) {
                String[] s = line.split(" ");

                // Gets file depending on format version.
                String file = (oldFormat ? s[0] : s[indexes[INDEX_IMAGE]]);

                ImagePlus imptmp = tr.openImage(directory, file);

                if (is == null) {
                    is = new ImageStack(imptmp.getWidth(), imptmp.getHeight());
                }
                is.addSlice(file, imptmp.getProcessor());
            }
        } catch (Exception e) {
            IJ.showStatus("");
            IJ.showMessage("Sel_Reader", "Sel_Reader : " + e);
            return null;
        }

        ImagePlus imp = new ImagePlus(fileName, is);

        return imp;
    }

    private static int getFieldIndex(String field) {
        if (field.toLowerCase().compareTo(IMAGE) == 0) {
            return INDEX_IMAGE;
        } else if (field.toLowerCase().compareTo(ENABLED) == 0) {
            return INDEX_ENABLED;
        }

        return -1;
    }

    /**
     *  Loads only the files list.
     *
     *@param  directory  Description of the Parameter
     *@param  fileName   Description of the Parameter
     *@return            Description of the Return Value
     */
    public static String[] loadFileNames(String directory, String fileName) {
        if ((fileName == null) || (fileName.isEmpty())) {
            return null;
        }

        if (!directory.endsWith(File.separator)) {
            directory += File.separator;
        }

        ArrayList<String> filesList = new ArrayList<String>();

        try {
            RandomAccessFile f = new RandomAccessFile(directory + fileName, "r");
            String line;
            boolean oldFormat;
            int indexes[] = null;

            // Skips header
            line = f.readLine();//; XMIPP_3 * column_format *
            if (line.toUpperCase().contains("XMIPP_3")) {
                oldFormat = false;

                f.readLine();// 2nd header line ";"

                // 3rd line contains parameters order.
                String metadata[] = f.readLine().split(" ");
                indexes = new int[metadata.length - 1];
                for (int i = 0; i < indexes.length; i++) {
                    indexes[i] = -1;
                }

                for (int i = 1; i < metadata.length; i++) { // 0 position is ";"
                    int index = getFieldIndex(metadata[i]);
                    if (index >= 0) {
                        indexes[index] = i - 1;
                    }
                }
            } else {  // Old format.
                oldFormat = true;
                f.seek(0);
            }

            while ((line = f.readLine()) != null) {
                String[] s = line.split(" ");

                // Gets file depending on format version.
                String file = (oldFormat ? s[0] : s[indexes[INDEX_IMAGE]]);

                if (!file.trim().isEmpty()) {
                    // Fixes file name by adding its parent.
                    if (!file.startsWith(File.separator)) {
                        file = directory /*+ File.separator */ + file.split(" ")[0];
                    }

                    filesList.add(file);
                }
            }
        } catch (Exception e) {
            IJ.showStatus("");
            IJ.showMessage("Sel_Reader", "Sel_Reader : " + e);
            return null;
        }

        return filesList.toArray(new String[filesList.size()]);
    }
}
