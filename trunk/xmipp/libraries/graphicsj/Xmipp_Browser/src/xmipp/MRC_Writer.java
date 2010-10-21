package xmipp;



import ij.plugin.filter.*;
import ij.*;
import ij.io.*;
import ij.process.*;
import java.io.*;

/**
 *  Description of the Class
 *
 *@author     cedric
 *@created    04/01/2005
 */
public class MRC_Writer implements PlugInFilter {
	ImagePlus img;
	String directory;
	String name;
	boolean littleEndian=true;

	/**
	 *  Main processing method for the SpiderReader_ object
	 *
	 *@param  arg  Description of the Parameter
	 *@param  imp  Description of the Parameter
	 *@return      Description of the Return Value
	 */
	public int setup(String arg, ImagePlus imp) {
		img = imp;
		if(!(arg==null||arg.equals(""))){
			File dest = new File(arg);
			directory = dest.getParent()+ System.getProperty("file.separator");
			name = dest.getName();
			System.out.println("saving as mrc "+directory+name);
		}else{
			name=null;
		}
		return DOES_32 + DOES_16 + DOES_8G + NO_CHANGES;
	}


	/**
	 *  Main processing method for the SpiderWriter_ object
	 *
	 *@param  ip  Description of the Parameter
	 */
	public void run(ImageProcessor ip) {
		if(name==null){
			SaveDialog sd = new SaveDialog("save as...", img.getTitle(), ".mrc");
			directory=sd.getDirectory();
			name=sd.getFileName();
		}
		if (name == null) {
			return;
		}
		IJ.showStatus("saving: " + directory + name);
		//ImageStack stack = myimp.getStack();
		//save(stack, directory + name);
		//SPIDER.saveImage(directory, name, myimp, true);


		byte[] header = createHeader();
		FileInfo fi = img.getFileInfo();
		fi.fileFormat = fi.RAW;
		System.out.println("saving...intel order=" + littleEndian);
		fi.intelByteOrder = littleEndian;
		fi.fileName = name;
		fi.directory = directory;
		ImageWriter file = new ImageWriter(fi);
		try {
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(directory +System.getProperty("file.separator")+ name)));
			for (byte aHeader : header) {
				out.writeByte(aHeader);
			}
			file.write(out);
			out.close();
		} catch (IOException ioe) {
			IJ.error("" + ioe);
		}
	}



	/**
	 *  Description of the Method
	 *
	 *@return    Description of the Return Value
	 */
	private byte[] createHeader() {
		ImageProcessor ip = img.getProcessor();
		int mode = 0;
		if (ip instanceof ByteProcessor) {
			mode = 0;
		} else if (ip instanceof ShortProcessor) {
			mode = 1;
		} else if (ip instanceof FloatProcessor) {
			mode = 2;
		}

		double sum = 0;
		int nbPix = 0;
		float min = Float.NaN;
		float max = Float.NaN;
		float value;
		int mini = 0;
		int maxi = 0;
		int avgi = 0;
		for (int z = 0; z < img.getNSlices(); z++) {
			ip = img.getImageStack().getProcessor(z + 1);
			for (int y = 0; y < img.getHeight(); y++) {
				for (int x = 0; x < img.getWidth(); x++) {
					value = ip.getPixelValue(x, y);
					sum += value;
					nbPix++;
					if (Float.isNaN(min)) {
						min = value;
					}
					if (Float.isNaN(max)) {
						max = value;
					}
					if (max < value) {
						max = value;
					} else if (min > value) {
						min = value;
					}
				}
			}
		}
		float avg = (float) (sum / nbPix);
		maxi = Float.floatToRawIntBits(max);
		mini = Float.floatToRawIntBits(min);
		avgi = Float.floatToRawIntBits(avg);

		byte[] header = new byte[1024];
		writeIntInHeader(img.getWidth(), 0 * 4, header, littleEndian);
		writeIntInHeader(img.getHeight(), 1 * 4, header, littleEndian);
		writeIntInHeader(img.getNSlices(), 2 * 4, header, littleEndian);
		writeIntInHeader(mode, 3 * 4, header, littleEndian);
		writeIntInHeader(1, 16 * 4, header, littleEndian);
		writeIntInHeader(2, 17 * 4, header, littleEndian);
		writeIntInHeader(3, 18 * 4, header, littleEndian);
		writeIntInHeader(mini, 19 * 4, header, littleEndian);
		writeIntInHeader(maxi, 20 * 4, header, littleEndian);
		writeIntInHeader(avgi, 21 * 4, header, littleEndian);
		return header;
	}


	/**
	 *  Description of the Method
	 *
	 *@param  value         Description of the Parameter
	 *@param  index         Description of the Parameter
	 *@param  h             Description of the Parameter
	 *@param  littleEndian  Description of the Parameter
	 */
	private static void writeIntInHeader(int value, int index, byte[] h, boolean littleEndian) {
		byte b1 = (byte) (value >> 24);
		byte b2 = (byte) ((value << 8) >> 24);
		byte b3 = (byte) ((value << 16) >> 24);
		byte b4 = (byte) ((value << 24) >> 24);
		if (littleEndian) {
			h[index] = b4;
			h[index + 1] = b3;
			h[index + 2] = b2;
			h[index + 3] = b1;
		} else {
			h[index] = b1;
			h[index + 1] = b2;
			h[index + 2] = b3;
			h[index + 3] = b4;
		}
	}

}
