import ij.plugin.*;
import ij.*;
import ij.io.*;
import java.io.*;
import java.nio.*;

/**
 *  Description of the Class
 *
 *@author     cedric
 *@created    04/01/2005
 */
public class Spider_Reader extends ImagePlus implements PlugIn {
	float[] header;
	boolean littleEndian = false;


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
			OpenDialog od = new OpenDialog("Load SPIDER File...", "");
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


	/**
	 *  Description of the Method
	 *
	 *@param  directory  Description of the Parameter
	 *@param  fileName   Description of the Parameter
	 *@return            Description of the Return Value
	 */
	public ImagePlus load(String directory, String fileName) {
		/*
		 *  throws IOException
		 */
		if ((fileName == null) || (fileName.equals(""))) {
			return null;
		}

		if (!directory.endsWith(File.separator)) {
			directory += File.separator;
		}

		IJ.showStatus("Loading SPIDER File: " + directory + fileName);

		// Try calling the parse routine
		try {
			parseSPIDER(directory, fileName);
		} catch (Exception e) {
			IJ.showStatus("parseSPIDER() error");
			IJ.showMessage("SPIDER_Reader", "" + e);
			return null;
		}

		// Set the file information
		FileInfo fi = new FileInfo();
		// Go and fetch the DM3 specific file Information
		try {
			fi = getSPIDERFileInfo(directory, fileName);
		}
		// This is in case of trouble parsing the tag table
		catch (Exception e) {
			IJ.showStatus("");
			IJ.showMessage("SPIDER_Reader", "gSPIDER:" + e);
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
	void parseSPIDER(String directory, String fileName) throws IOException {
		// This reads through the MRC file, extracting useful tags
		// which allow one to determine the data offset etc.

		RandomAccessFile stream = new RandomAccessFile(directory + fileName, "r");
		//IJ.write("test");
		stream.seek(4 * 4);
		float iform = stream.readFloat();
		//System.out.println("iform" + iform);
		//stream.seek(0);
		if (iform != 1 && iform != 3 && iform != -11 && iform != -12 && iform != -21 && iform != -22) {
			//IJ.write("little endian byte order");
			littleEndian = true;
			//System.out.println("little endian:" + littleEndian);
			iform = changeByteOrder(iform);
			//System.out.println("iform" + iform);
			if (iform != 1 && iform != 3 && iform != -11 && iform != -12 && iform != -21 && iform != -22) {
				throw new IOException("cannot determine byte order in SPIDER file.");
			}
			//stream.seek(0);
		}

		stream.seek(21 * 4);
		float size = stream.readFloat();
		if (littleEndian) {
			size = changeByteOrder(size);
		}
		int hsize = (int) size;
		stream.seek(0);
		header = new float[hsize / 4];
		byte[] tmp = new byte[hsize];
		stream.readFully(tmp, 0, hsize);
		//System.out.println("all header taken from file");
		for (int i = 0; i < hsize / 4; i++) {
			int orig;
			if (littleEndian) {
				orig = ((((tmp[i * 4 + 3] & 0xff) << 24) | ((tmp[i * 4 + 2] & 0xff) << 16) | ((tmp[i * 4 + 1] & 0xff) << 8) | (tmp[i * 4] & 0xff)));
			} else {
				orig = ((((tmp[i * 4] & 0xff) << 24) | ((tmp[i * 4 + 1] & 0xff) << 16) | ((tmp[i * 4 + 2] & 0xff) << 8) | (tmp[i * 4 + 3] & 0xff)));
			}
			header[i] = Float.intBitsToFloat(orig);
		}
		// Close the input stream
		stream.close();
	}


	/**
	 *  Gets the mRCFileInfo attribute of the MRC_Reader object
	 *
	 *@param  directory        Description of the Parameter
	 *@param  fileName         Description of the Parameter
	 *@return                  The mRCFileInfo value
	 *@exception  IOException  Description of the Exception
	 */
	FileInfo getSPIDERFileInfo(String directory, String fileName) throws IOException {
		// this gets the basic file information using the contents of the tag
		// tables created by parseMRC()
		FileInfo fi = new FileInfo();
		fi.fileFormat = fi.RAW;
		fi.fileName = fileName;
		fi.directory = directory;
		// Originally forgot to do this - tells ImageJ what endian form the actual
		// image data is in
		if (littleEndian) {
			fi.intelByteOrder = littleEndian;
		}
		fi.fileType = FileInfo.GRAY32_FLOAT;
		fi.width = (int) header[11];
		fi.height = (int) header[1];
		fi.nImages = (int) header[0];
		fi.offset = (int) header[21];
		return fi;
	}


	/**
	 *  Description of the Method
	 *
	 *@param  value  Description of the Parameter
	 *@return        Description of the Return Value
	 */
	float changeByteOrder(float value) {
		ByteBuffer bb = ByteBuffer.allocateDirect(4);
		bb.putFloat(0, value);
		bb.order(ByteOrder.LITTLE_ENDIAN);
		return bb.getFloat(0);
	}

}
