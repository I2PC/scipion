import ij.plugin.*;
import ij.*;
import ij.io.*;
import ij.process.*;
import java.io.*;

/**
 *  Reading MRC files modified the 18/10/2007 adding substraction of 32768 to
 *  the pixels value of FEI files - FEI don't use correct version of MRC format!
 *
 *@author     cedric
 *@created    04/01/2005
 */
public class MRC_Reader extends ImagePlus implements PlugIn {
	byte[] header;
	boolean FEI = false;


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
			OpenDialog od = new OpenDialog("Load MRC File...", arg);
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
		if (FEI) {
			ImageStack stack = imp.getStack();
			for (int i = 1; i <= stack.getSize(); i++) {
				stack.getProcessor(i).add(-32768);
			}
			imp.getProcessor().resetMinAndMax();
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
			parseMRC(directory, fileName);
		} catch (Exception e) {
			IJ.showStatus("parseMRC() error");
			IJ.showMessage("DM3_Reader", "" + e);
			return null;
		}

		// Set the file information
		FileInfo fi = new FileInfo();
		// Go and fetch the DM3 specific file Information
		try {
			fi = getMRCFileInfo(directory, fileName);
		}
		// This is in case of trouble parsing the tag table
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
	 *@exception  IOException  Description of the Exception
	 */
	void parseMRC(String directory, String fileName) throws IOException {
		// This reads through the MRC file, extracting useful tags
		// which allow one to determine the data offset etc.

		// alternative way to read from file - allows seeks!
		// and therefore keeps track of position (use long getFilePointer())
		// also has DataInput Interface allowing
		// reading of specific types
		RandomAccessFile f = new RandomAccessFile(directory + fileName, "r");

		header = new byte[1024];
		for (int i = 0; i < 1024; i++) {
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
	FileInfo getMRCFileInfo(String directory, String fileName) throws IOException {
		// this gets the basic file information using the contents of the tag
		// tables created by parseMRC()
		FileInfo fi = new FileInfo();
		fi.fileFormat = fi.RAW;
		fi.fileName = fileName;
		fi.directory = directory;
		// Originally forgot to do this - tells ImageJ what endian form the actual
		// image data is in
		boolean littleEndian = false;
		int mode = readIntInHeader(3 * 4, header, littleEndian);
		if (mode == 0) {
			fi.fileType = FileInfo.GRAY8;
			fi.width = readIntInHeader(0, header, littleEndian);
			fi.height = readIntInHeader(1 * 4, header, littleEndian);
			fi.nImages = readIntInHeader(2 * 4, header, littleEndian);
			if (fi.width < 0 || fi.height < 0 || fi.nImages < 0) {
				littleEndian = true;
			}
		} else if (mode == 1) {
			fi.fileType = FileInfo.GRAY16_UNSIGNED;
		} else if (mode == 2) {
			fi.fileType = FileInfo.GRAY32_FLOAT;
		} else {
			//it is not a big endian data
			littleEndian = true;
			mode = readIntInHeader(3 * 4, header, littleEndian);
			if (mode == 1) {
				fi.fileType = FileInfo.GRAY16_UNSIGNED;
			} else if (mode == 2) {
				fi.fileType = FileInfo.GRAY32_FLOAT;
			} else {
				throw new IOException("Unimplemented ImageData mode=" + mode + " in MRC file.");
			}
		}
		//littleEndian = true;
		if (littleEndian) {
			fi.intelByteOrder = littleEndian;
		}
		fi.width = readIntInHeader(0, header, littleEndian);
		fi.height = readIntInHeader(1 * 4, header, littleEndian);
		fi.nImages = readIntInHeader(2 * 4, header, littleEndian);
		int extrabytes = readIntInHeader(23 * 4, header, littleEndian);
		fi.offset = 1024 + extrabytes;
		//IJ.write("" + extrabytes);
		if (extrabytes == 131072) {
			IJ.write("FEI");
			//fi.fileType = FileInfo.GRAY16_SIGNED;
			FEI = true;
		}

		IJ.write("mode=" + mode + " littleEndian=" + littleEndian);
		IJ.write("x=" + fi.width + " y=" + fi.height + " z=" + fi.nImages);
		IJ.write("extrabytes=" + extrabytes);

		return fi;
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
