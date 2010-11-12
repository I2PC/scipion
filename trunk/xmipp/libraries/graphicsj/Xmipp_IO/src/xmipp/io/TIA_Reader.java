package xmipp.io;

/*
 * Modified to extend ImagePlus and to allow total integration.
 * @author Juanjo Vega
 */
/*

This plugin tries to import the *.ser-files (STEM, CCD) written by the TIA Program (Emispec) during the data aquisition. I used the on Dr. Chris Boothroyd website described file structure (http://www.microscopy.cen.dtu.dk/~cbb/info/TIAformat/index.html). I want to thank him for the given information on his website.

@author Steffen Schmidt; steffen.schmidt at cup.uni-muenchen.de

In this plugin the ledatastream class of the Canadian Mind Products company is used (please see the additinal license rule).

The drag and drop support in ImageJ can be enabled by editing/compiling the HandleExtraFileTypes.java-file in the Input-Output folder. I added the line: if (name.endsWith(".ser")) {return tryPlugIn("TIA_Reader", path); in the tryOpen-class and recompiled it.

 */
import com.mindprod.ledatastream.LEDataInputStream;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.PlotWindow;
import ij.io.FileInfo;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import xmipp.io.ij.FileOpener;

public class TIA_Reader extends ImagePlus implements PlugIn {

    int byteoffset;

    public void run(String arg) {
        String path = getPath(arg);
        if (null == path) {
            return;
        }
        if (!parseTIA(path)) {
            return;
        }
    }

    private String getPath(String arg) {
        if (null != arg) {
            if (0 == arg.indexOf("http://")
                    || new File(arg).exists()) {
                return arg;
            }
        }
        OpenDialog od = new OpenDialog(".ser-Reader...", null);
        String DIRECTORY = od.getDirectory();
        if (null == DIRECTORY) {
            return null; // dialog was canceled
        }
        DIRECTORY = DIRECTORY.replace('\\', '/'); // Windows-unix convention
        if (!DIRECTORY.endsWith("/")) {
            DIRECTORY += "/";
        }
        String FILENAME = od.getFileName();
        if (null == FILENAME) {
            return null; //dialog was canceled
        }
        return DIRECTORY + FILENAME;
    }

    private InputStream open(String path) throws Exception {
        if (0 == path.indexOf("http://")) {
            return new java.net.URL(path).openStream();
        }
        return new FileInputStream(path);
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

        IJ.showStatus("Loading SPIDER File: " + directory + fileName);

        // Try calling the parse routine
        try {
            parseTIA(directory + fileName);
        } catch (Exception e) {
            IJ.showStatus("parseSPIDER() error");
            IJ.showMessage("SPIDER_Reader", "" + e);
            return null;
        }

        // Set the file information
        FileInfo fi = new FileInfo();
        // Go and fetch the DM3 specific file Information
        try {
            fi = getTIAFileInfo(directory + fileName);
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

    public boolean parseTIA(String path) {

        //variables
        int FILE_OFFSET; //number of bytes for jumping to the data
        int NUMBER_IMAGES; //number data sets
        int OFFSET_ARRAY_OFFSET; //Offset to the data array offset
        int[] DATA_OFFSET; //field of data offset array values
        int DATA_TYPE_ID; //type of stored data 0x4120->1D; 0x4122->2D

        //reading the header and the data offset array
        try {
            InputStream is = open(path);
            LEDataInputStream data = new LEDataInputStream(is);
            if (data.readShort() != 0x4949) {
                IJ.error("no little-endian byte ordering");
                return false; //ByteOrder 18761=0x4949H indicates little-endian byte Ordering
            }
            data.readShort(); //SeriesID
            data.readShort(); //SeriesVersion
            DATA_TYPE_ID = data.readInt(); //DataTypeID
            data.readInt(); //TagTypeID
            data.readInt(); //TotalNumberElements
            NUMBER_IMAGES = data.readInt(); //ValidNumberElements
            OFFSET_ARRAY_OFFSET = data.readInt(); //OffsetArrayOffset
            data.readInt(); //NumberDimension
            DATA_OFFSET = new int[NUMBER_IMAGES]; //configure the size of the data offset array
            data.skipBytes(OFFSET_ARRAY_OFFSET - 30); //Data offset array - header format
            int count = 0;
            while (count < NUMBER_IMAGES) {
                DATA_OFFSET[count] = data.readInt();
                count++;
            } //reading data offset values
            data.close();
        } catch (Exception e) {
            IJ.error("Error opening file");
            return false;
        }

        //opening of the different data elements
        try {
            int count = 0;
            while (count < NUMBER_IMAGES) {
                if (DATA_TYPE_ID == 0x4122) {
                    OpenImage(path, DATA_OFFSET[count]);
                } //DataTypeID = 0x4122 indicates an image
                else if (DATA_TYPE_ID == 0x4120) {
                    OpenSpectra(path, DATA_OFFSET[count]);
                } //DataTypeID = 0x4120 indicates an spectra
                else if (check_data_element(path, DATA_OFFSET[count])) //guessing of the DataType
                {
                    OpenImage(path, DATA_OFFSET[count]);
                } else {
                    OpenSpectra(path, DATA_OFFSET[count]);
                }
                count++;
            }
        } catch (Exception e) {
            IJ.error("Error opening Data series");
            return false;
        }
        return true;
    }

    private boolean check_data_element(String path, int byteoffset) throws Exception {

        //variables
        double PIXEL_WIDTH; //CalibrationDeltaX
        double PIXEL_HEIGHT; //CalibrationDeltaY
        short DATA_TYPE; //DataType

        // reading the header of the data elements
        InputStream is = open(path);
        LEDataInputStream data = new LEDataInputStream(is);
        data.skipBytes(byteoffset); //jumping to the data element field
        data.readDouble(); //CalibrationOffsetX
        PIXEL_WIDTH = data.readDouble(); //CalibrationDeltaX
        data.readInt();  //CalibrationElementX
        DATA_TYPE = data.readShort(); //may be DataType
        data.skipBytes(6); //jumping to the CalibrationDeltaY
        PIXEL_HEIGHT = data.readDouble(); //CalibrationDeltaY
        data.close(); //close file

        //guessing of the DataType; true indicates 2D - false indicates 1D
        if (DATA_TYPE == 0) {
            return true;
        } else if (PIXEL_WIDTH == PIXEL_HEIGHT) {
            return true;
        } else {
            return false;
        }
    }

    public FileInfo getTIAFileInfo(String path) throws Exception {

        //variables
        double PIXEL_WIDTH; //CalibrationDeltaX
        double PIXEL_HEIGHT; //CalibrationDeltaY
        short DATA_TYPE; //DataType
        int IMAGE_WIDTH; //ArraySizeX
        int IMAGE_HEIGHT; //ArraySizeY

        // reading calibration values
        InputStream is = open(path);
        LEDataInputStream data = new LEDataInputStream(is);
        data.skipBytes(byteoffset); //jumping to the 2D-data element field
        data.readDouble(); //CalibrationOffsetX
        PIXEL_WIDTH = data.readDouble(); //CalibrationDeltaX
        data.readInt();  //CalibrationElementX
        data.readDouble(); //CalibrationOffsetY
        PIXEL_HEIGHT = data.readDouble(); //CalibrationDeltaY
        data.readInt(); //CalibrationElementY
        DATA_TYPE = data.readShort(); //DataType
        IMAGE_WIDTH = data.readInt(); //ArraySizeX
        IMAGE_HEIGHT = data.readInt(); //ArraySizeY
        data.close();

//					IJ.log ("Dispersion-x:  "+(PIXEL_WIDTH)+" / Dispersion-y: "+(PIXEL_HEIGHT));
//					IJ.log ("Array-x: "+(IMAGE_WIDTH)+" / Array-y: "+(IMAGE_HEIGHT));

        //opening of the image
        FileInfo fi = new FileInfo();
        fi.fileFormat = fi.RAW;
        fi.intelByteOrder = true;  //little-endian byte ordering
        int islash = path.lastIndexOf('/');
        if (0 == path.indexOf("http://")) {
            fi.url = path;
        } else {
            fi.directory = path.substring(0, islash + 1);
        }
        fi.fileName = path.substring(islash + 1);
        fi.width = IMAGE_WIDTH;
        fi.height = IMAGE_HEIGHT;
        if ((PIXEL_WIDTH * IMAGE_WIDTH) < 1E-5) {
            fi.pixelWidth = PIXEL_WIDTH / 1E-9;
            fi.pixelHeight = PIXEL_HEIGHT / 1E-9;
            fi.unit = "nm";
        } else if ((PIXEL_WIDTH * IMAGE_WIDTH) < 1E-2) {
            fi.pixelWidth = PIXEL_WIDTH / 1E-6;
            fi.pixelHeight = PIXEL_HEIGHT / 1E-6;
            fi.unit = "microns";
        } else if ((PIXEL_WIDTH * IMAGE_WIDTH) < 1E1) {
            fi.pixelWidth = PIXEL_WIDTH / 1E-3;
            fi.pixelHeight = PIXEL_HEIGHT / 1E-3;
            fi.unit = "mm";
        } else {
            fi.pixelWidth = PIXEL_WIDTH;
            fi.pixelHeight = PIXEL_HEIGHT;
            fi.unit = "m";
        }
        fi.nImages = 1;
        switch (DATA_TYPE) {
            case 1:
                fi.fileType = FileInfo.GRAY8;
                break;
            case 2:
                fi.fileType = FileInfo.GRAY16_UNSIGNED;
                break;
            case 3:
                fi.fileType = FileInfo.GRAY32_UNSIGNED;
                break;
            case 4:
                fi.fileType = FileInfo.GRAY8;
                break;
            case 5:
                fi.fileType = FileInfo.GRAY16_SIGNED;
                break;
            case 6:
                fi.fileType = FileInfo.GRAY32_INT;
                break;
            case 7:
                fi.fileType = FileInfo.GRAY32_FLOAT;
                break;
            case 8:
                fi.fileType = FileInfo.GRAY64_FLOAT;
                break;
        }
        fi.offset = byteoffset + 50;
        fi.whiteIsZero = false;

        return fi;
    }

    private void OpenImage(String path, int byteoffset) throws Exception {
        this.byteoffset = byteoffset;
        FileInfo fi = getTIAFileInfo(path);

        FileOpener fo = new FileOpener(fi);
        ImagePlus imp = fo.open(false);
        IJ.run(imp, "Flip Vertically", "");
        Object obinfo = imp.getProperty("Info");
        if (null != obinfo) {
            imp.setProperty("Info", obinfo);
        }

        // Sets "this" data with loaded stuff.
        if (imp != null) {
            String fileName = (new File(path)).getName();

            // Attach the Image Processor
            if (imp.getNSlices() > 1) {
                setStack(fileName, imp.getStack());
            } else {
                setProcessor(fileName, imp.getProcessor());
            }
            // Copy the scale info over
            copyScale(imp);
        }
    }

    private void OpenSpectra(String path, int byteoffset) throws Exception {
        //variables
        double PIXEL_WIDTH; //CalibrationDelta
        double CALIBRATION_OFFSET; //CalibrationOffset
        int CALIBRATION_ELEMENT; //CalibrationElement
        short DATA_TYPE; //DataType
        int IMAGE_WIDTH; //ArrayLength
        float[] x; //x-value
        float[] y; //y-value

        // reading the calibration data
        InputStream is = open(path);
        LEDataInputStream data = new LEDataInputStream(is);
        data.skipBytes(byteoffset); //jumping to the data element
        CALIBRATION_OFFSET = data.readDouble(); //CalibrationOffset
        PIXEL_WIDTH = data.readDouble(); //CalibrationDelta
        CALIBRATION_ELEMENT = data.readInt();  //CalibrationElement
        DATA_TYPE = data.readShort(); //DataType
        IMAGE_WIDTH = data.readInt(); //ArrayLength
        x = new float[IMAGE_WIDTH]; //x-values
        y = new float[IMAGE_WIDTH]; //y-values

        //opening of the spectra
        int count = 0;
        while (count < IMAGE_WIDTH) {
            x[count] = 0.001f * ((float) CALIBRATION_OFFSET - ((float) PIXEL_WIDTH * (float) CALIBRATION_ELEMENT) + ((float) count * (float) PIXEL_WIDTH)); //setting in x in keV
            switch (DATA_TYPE) {
                case 1:
                    return;
                case 2:
                    y[count] = data.readShort();
                    break;
                case 3:
                    y[count] = data.readInt();
                    break;
                case 4:
                    return;
                case 5:
                    y[count] = data.readShort();
                    break;
                case 6:
                    y[count] = data.readInt();
                    break;
                case 7:
                    y[count] = data.readFloat();
                    break;
                case 8:
                    return;
                case 9:
                    return;
                case 10:
                    return;
            }
            count++;
        }
        data.close();
//					IJ.log ("Calibration-Delta-x: " + (PIXEL_WIDTH) + " / Length of data field: " + (IMAGE_WIDTH));
        //plotting of the data
        int islash = path.lastIndexOf('/');
        PlotWindow plot = new PlotWindow(path.substring(islash + 1), "keV", "counts", x, y);
        plot.setLineWidth(1);
        plot.draw();
    }
}


