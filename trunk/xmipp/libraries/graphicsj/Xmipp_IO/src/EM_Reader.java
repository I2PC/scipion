

import ij.io.OpenDialog;
import ij.io.FileInfo;
import xmipp.ij.io.FileOpener;
import ij.ImagePlus;
import ij.IJ;
import ij.plugin.PlugIn;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;

/**
 * Created by IntelliJ IDEA.
 *
 * @author MESSAOUDI Cedric
 *         Date: 26 nov. 2009
 *         Time: 10:16:12
 */
public class EM_Reader extends ImagePlus implements PlugIn {

    byte[] header;
    boolean littleEndian = false;

    /**
     *  Main processing method for the SpiderReader_ object
     *
     *@param  arg  Description of the Parameter
     */
    public void run(String arg) {
        String directory = "";
        String fileName = arg;
        if ((arg == null) || (arg.isEmpty())) {
            // Choose a file since none specified
            OpenDialog od = new OpenDialog("Load EM File...", arg);
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
        setProperties(header);

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
    public ImagePlus load(String directory, String fileName) {
        /*
        throws IOException
         */
        if ((fileName == null) || (fileName.equals(""))) {
            return null;
        }

        if (!directory.endsWith(File.separator)) {
            directory += File.separator;
        }

        IJ.showStatus("Loading MRC File: " + directory + fileName);

        // Try calling the parse routine
        try {
            parseEM(directory, fileName);
        } catch (Exception e) {
            IJ.showStatus("parseMRC() error");
            IJ.showMessage("DM3_Reader", "" + e);
            return null;
        }

        // Set the file information
        FileInfo fi = new FileInfo();
        // Go and fetch the DM3 specific file Information
        try {
            fi = getEMFileInfo(directory, fileName);
        } // This is in case of trouble parsing the tag table
        catch (Exception e) {
            IJ.showStatus("");
            IJ.showMessage("MRC_Reader", "gMRC:" + e);
            return null;
        }

        // Open the image!
        FileOpener fo = new FileOpener(fi);
        return fo.open(false);
    }

    /**
     *  Description of the Method
     *
     *@param  directory        Description of the Parameter
     *@param  fileName         Description of the Parameter
     *@exception java.io.IOException  Description of the Exception
     */
    void parseEM(String directory, String fileName) throws IOException {
        // This reads through the EM file, extracting useful tags
        // which allow one to determine the data offset etc.

        // alternative way to read from file - allows seeks!
        // and therefore keeps track of position (use long getFilePointer())
        // also has DataInput Interface allowing
        // reading of specific types
        RandomAccessFile f = new RandomAccessFile(directory + fileName, "r");

        header = new byte[512];
        for (int i = 0; i < 512; i++) {
            header[i] = f.readByte();
        }

        // Close the input stream
        f.close();
    }

    /**
     *  Gets the mRCFileInfo attribute of the MRC_Reader object
     *
     *@param  directory        Description of the Parameter
     *@param  fileName         Description of the Parameter
     *@return                  The mRCFileInfo value
     *@exception  IOException  Description of the Exception
     */
    FileInfo getEMFileInfo(String directory, String fileName) throws IOException {
        // this gets the basic file information using the contents of the tag
        // tables created by parseMRC()
        FileInfo fi = new FileInfo();
        fi.fileFormat = fi.RAW;
        fi.fileName = fileName;
        fi.directory = directory;
        // Originally forgot to do this - tells ImageJ what endian form the actual
        // image data is in
        littleEndian = header[0] > 4;
        System.out.println("little endian=" + littleEndian);
        int mode = header[3];
        String type = "";
        switch (mode) {
            case 1:
                fi.fileType = FileInfo.GRAY8;
                type += "byte";
                break;
            case 2:
                fi.fileType = FileInfo.GRAY16_UNSIGNED;
                type += "short unsigned";
                break;
            case 4:
                fi.fileType = FileInfo.GRAY32_INT;
                type += "int signed";
                break;
            case 5:
                fi.fileType = FileInfo.GRAY32_FLOAT;
                type += "float";
                break;
            case 9:
            default:
                fi.fileType = FileInfo.GRAY64_FLOAT;
                type += "double";
        }
        fi.intelByteOrder = littleEndian;

        fi.width = readIntInHeader(4, header, littleEndian);
        fi.height = readIntInHeader(8, header, littleEndian);
        fi.nImages = readIntInHeader(12, header, littleEndian);
        fi.offset = 512;
        System.out.println("type=" + type);
        System.out.println("width=" + fi.width);
        System.out.println("height=" + fi.height);
        System.out.println("depth=" + fi.nImages);
        System.out.println("offset=" + fi.offset);




        return fi;
    }

    /**
     * set the properties of image found in the header
     * <ul>
     * <li> key:unit:significance</li>
     * <li>U:Volt:accelerating volatage</li>
     * <li>COE:???m:Cs of objective lense</li>
     * <li>APE:mrad:aperture</li>
     * <li>VE:x:end magnification</li>
     * <li>VN:-:postmagnification of CCD</li>
     * <li>ET:s:exposure time in seconds</li>
     * <li>XG:-:pixel size in object-plane</li>
     * <li>DG:nm:EM-Code:EM420=1;CM12=2;CM200=3;CM120/Biofileter=4;CM300=5;Polara=6;Titan=7;extern=0</li>
     * <li>APD:nm:photometer aperture</li>
     * <li>L:nm:phys_pixel_size*nr_of_pixels</li>
     * <li>DF:Angstr:defocus, underfocus is neg</li>
     * <li>FA:Angstr:astigmatism</li>
     * <li>PHI:deg/1000:angle of astigmatism</li>
     * <li>DR:Angstr:drift</li>
     * <li>DELT:deg/1000:direction of drift</li>
     * <li>DDF:Angstr:focusincr for focus-series</li>
     * <li>KW:deg/1000:tiltangle</li>
     * <li>KR:deg/1000:axis perpend to tiltaxis</li>
     * </ul>
     *to retrieve these informations use the method getProperty with the corresponding key
     * @param header header of the em file
     * @see ImagePlus.getProperty(String key)
     */
    public void setProperties(byte[] header) {
        setProperty("U", readIntInHeader(96, header, littleEndian) / 1000000.0);
        setProperty("COE", readIntInHeader(100, header, littleEndian) / 1000.0);
        setProperty("APE", readIntInHeader(104, header, littleEndian) / 1000.0);
        setProperty("VE", readIntInHeader(108, header, littleEndian));
        setProperty("VN", readIntInHeader(112, header, littleEndian) / 1000.0);
        setProperty("ET", readIntInHeader(116, header, littleEndian) / 1000.0);
        setProperty("XG", readIntInHeader(120, header, littleEndian));
        setProperty("DG", readIntInHeader(124, header, littleEndian) / 1000.0);
        setProperty("APD", readIntInHeader(128, header, littleEndian) / 1000.0);
        setProperty("L", readIntInHeader(132, header, littleEndian) / 1000.0);
        setProperty("DF", readIntInHeader(136, header, littleEndian));
        setProperty("FA", readIntInHeader(140, header, littleEndian));
        setProperty("PHI", readIntInHeader(144, header, littleEndian) / 1000.0);
        setProperty("DR", readIntInHeader(148, header, littleEndian));
        setProperty("DELT", readIntInHeader(152, header, littleEndian) / 1000.0);
        setProperty("DDF", readIntInHeader(156, header, littleEndian));
        setProperty("KW", readIntInHeader(168, header, littleEndian) / 1000.0);
        setProperty("KR", readIntInHeader(172, header, littleEndian) / 1000.0);

        System.out.println("property added voltage kV:" + getProperty("U"));
        System.out.println("property added end magnification:" + getProperty("VE"));
        System.out.println("property added post magnification CCD:" + getProperty("VN"));
        System.out.println("property added exp time (s):" + getProperty("ET"));
        System.out.println("property added pixelsize:" + getProperty("XG"));
        System.out.println("property added physpixsize*nr of pixels:" + getProperty("L"));
    }

    /**
     *  Description of the Method
     *
     *@param  index         Description of the Parameter
     *@param  h             Description of the Parameter
     *@param  littleEndian  Description of the Parameter
     *@return               Description of the Return Value
     */
    private static int readIntInHeader(int index, byte[] h, boolean littleEndian) {
        byte b1 = h[index];
        byte b2 = h[index + 1];
        byte b3 = h[index + 2];
        byte b4 = h[index + 3];
        if (littleEndian) {
            return ((((b4 & 0xff) << 24) | ((b3 & 0xff) << 16) | ((b2 & 0xff) << 8) | (b1 & 0xff)));
        }
        return ((((b1 & 0xff) << 24) | ((b2 & 0xff) << 16) | ((b3 & 0xff) << 8) | (b4 & 0xff)));
    }
}
