package xmipp.io;

/**
 * Modified by
 * @author Juanjo Vega
 */
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileInfo;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import xmipp.io.ij.FileOpener;

public class OpenSPE_ extends ImagePlus implements PlugIn {

    public void run(String arg) {
        String directory = "";
        String fileName = arg;
        if ((arg == null) || (arg.compareTo("") == 0)) {
            // Choose a file since none specified
            OpenDialog od = new OpenDialog("Load SPE File...", "");
            fileName = od.getFileName();
            if (fileName == null) {
                return;
            }
            directory = od.getDirectory();
        } else {
            // we were sent a filename to openThumbnail
            File dest = new File(arg);
            directory = dest.getParent();
            fileName = dest.getName();
        }
        //System.out.println("" + directory + "       " + fileName);
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

    public static ImagePlus load(String directory, String file) {
        /*        File f = new File(directory, file);*/
        try {
            /*            FileInputStream in = new FileInputStream(f);
            byte[] h = new byte[SpeHeader.headerSize];
            in.read(h, 0, h.length);
            SpeHeader header = new SpeHeader(h);
            int speType = header.getDatatype();
            int fiType = fileInfoType(speType);
            if (fiType < 0) {
            IJ.showMessage("Open SPE...", "Invalid data type.");
            return null;
            }*/

            //if (speType == header.UNINT) {
            //	boolean convert = IJ.showMessageWithCancel("Open SPE...",
            //		"Convert UNSIGNED 16-bit integer \n to SIGNED 16-bit integer ?");
            //	if(convert)
            //	fiType = FileInfo.GRAY16_SIGNED;
            //}
//
//            FileInfo fi = new FileInfo();
//            fi.directory = directory;
//            fi.fileFormat = fi.RAW;
//            fi.fileName = file;
//            fi.fileType = fiType;
//            fi.gapBetweenImages = 0;
//            fi.height = header.getHeight();
//            fi.intelByteOrder = true;
//            fi.nImages = header.getStackSize();
//            fi.offset = SpeHeader.RAW_OFFSET;
//            fi.width = header.getWidth();
            FileInfo fi = getSPEFileInfo(directory, file);
            FileOpener fo = new FileOpener(fi);
            ImagePlus imp = fo.open(false);
            IJ.showStatus("");
            return imp;
        } catch (IOException e) {
            IJ.error("An error occured reading the file.\n \n" + e);
            IJ.showStatus("");
            return null;
        }
    }

    /**
     * @author Juanjo Vega
     *  Description of the Method
     *
     *@param  directory  Description of the Parameter
     *@param  fileName   Description of the Parameter
     *@return            Description of the Return Value
     */
    public ImagePlus loadThumbnail(String directory, String fileName, int w, int h) {
        return loadThumbnail(directory, fileName, w, h, FileOpener.MID_SLICE);
    }

    public ImagePlus loadThumbnail(String directory, String fileName, int w, int h, int slice) {
        if ((fileName == null) || (fileName.equals(""))) {
            return null;
        }

        if (!directory.endsWith(File.separator)) {
            directory += File.separator;
        }

        IJ.showStatus("Loading SPE File: " + directory + fileName);

        // Try calling the parse routine
        try {
            parseSPE(directory, fileName);
        } catch (Exception e) {
            IJ.showStatus("parseSPE() error");
            IJ.showMessage("SPE_Reader", "" + e);
            return null;
        }

        // Set the file information
        FileInfo fi = new FileInfo();
        // Go and fetch the SPE specific file Information
        try {
            fi = getSPEFileInfo(directory, fileName);
        } // This is in case of trouble parsing the tag table
        catch (Exception e) {
            IJ.showStatus("");
            IJ.showMessage("SPIDER_Reader", "gSPIDER:" + e);
            return null;
        }

        // Open the image!
        FileOpener fo = new FileOpener(fi);

        // Get mid slice if required.
        if (slice == FileOpener.MID_SLICE) {
            slice = (int) Math.ceil(Math.ceil((double) (fi.nImages / 2)));
        }

        // If thumbnail is bigger than the original image, loads it with the original size and zooms it.
        int th_size = w * h;
        int ip_size = fi.width * fi.height;

        ImagePlus thumbnail;
        if (th_size < ip_size) {
            thumbnail = fo.openThumbnail(false, slice, w, h);
        } else {
            thumbnail = fo.openThumbnail(false, slice, fi.width, fi.height);
            thumbnail.setProcessor(thumbnail.getTitle(), thumbnail.getProcessor().resize(w, h));
        }

        return thumbnail;
    }

    public static FileInfo getSPEFileInfo(String directory, String fileName) throws IOException {
        File f = new File(directory, fileName);

        FileInputStream in = new FileInputStream(f);
        byte[] h = new byte[SpeHeader.headerSize];
        in.read(h, 0, h.length);
        SpeHeader header = new SpeHeader(h);
        int speType = header.getDatatype();
        int fiType = fileInfoType(speType);
        if (fiType < 0) {
            IJ.showMessage("Open SPE...", "Invalid data type.");
            return null;
        }

        FileInfo fi = new FileInfo();
        fi.directory = directory;
        fi.fileFormat = fi.RAW;
        fi.fileName = fileName;
        fi.fileType = fileInfoType(speType);
        fi.gapBetweenImages = 0;
        fi.height = header.getHeight();
        fi.intelByteOrder = true;
        fi.nImages = header.getStackSize();
        fi.offset = SpeHeader.RAW_OFFSET;
        fi.width = header.getWidth();

        return fi;
    }

    public void parseSPE(String directory, String fileName) throws IOException {
        if (!directory.endsWith(File.separator)) {
            directory += File.separator;
        }
    }

    public static int fileInfoType(int speType) {
        switch (speType) {
            case SpeHeader.FLOAT:
                return FileInfo.GRAY32_FLOAT;
            case SpeHeader.LONG:
                return FileInfo.GRAY32_INT;
            case SpeHeader.INT:
                return FileInfo.GRAY16_SIGNED;
            case SpeHeader.UNINT:
                return FileInfo.GRAY16_UNSIGNED;
            default:
                return -1;
        }
    }
}
