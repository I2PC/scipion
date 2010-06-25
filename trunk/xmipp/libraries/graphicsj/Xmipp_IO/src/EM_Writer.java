import ij.plugin.filter.PlugInFilter;
import ij.ImagePlus;
import ij.IJ;
import ij.io.SaveDialog;
import ij.io.FileInfo;
import ij.io.ImageWriter;
import ij.process.ImageProcessor;
import ij.process.ByteProcessor;
import ij.process.ShortProcessor;
import ij.process.FloatProcessor;

import java.io.*;

/**
 * Created by IntelliJ IDEA.
 *
 * @author MESSAOUDI Cèdric
 *         Date: 26 nov. 2009
 *         Time: 12:15:27
 */
public class EM_Writer implements PlugInFilter {
	ImagePlus img;
	String directory;
	String name;
	boolean littleEndian=IJ.isWindows()||IJ.isLinux()||IJ.isMacOSX()||IJ.isVista();

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
			System.out.println("saving as em "+directory+name);
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
			SaveDialog sd = new SaveDialog("save as...", img.getTitle(), ".em");
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
		byte mode = 0;
		if (ip instanceof ByteProcessor) {
			mode = 1;
		} else if (ip instanceof ShortProcessor) {
			mode = 2;
		} else if (ip instanceof FloatProcessor) {
			mode = 5;
		}



		byte[] header = new byte[512];
		if (IJ.isWindows()||IJ.isVista()||IJ.isLinux()) header[0]=6;
		else if(IJ.isMacOSX()) header[0]=5;
		else if(IJ.isMacintosh()&&!IJ.isMacOSX()) header[0]=0;
		else header[0]=3;
		header[2]=0;
		header[3]=mode;
		writeIntInHeader(img.getWidth(), 4, header, littleEndian);
		writeIntInHeader(img.getHeight(), 8, header, littleEndian);
		writeIntInHeader(img.getNSlices(),12, header, littleEndian);

		
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
