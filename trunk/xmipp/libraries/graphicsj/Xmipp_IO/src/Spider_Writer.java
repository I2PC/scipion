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
public class Spider_Writer implements PlugInFilter {
	ImagePlus myimp;
	String directory;
	String name;

	/**
	 *  Main processing method for the SpiderReader_ object
	 *
	 *@param  arg  Description of the Parameter
	 *@param  imp  Description of the Parameter
	 *@return      Description of the Return Value
	 */
	public int setup(String arg, ImagePlus imp) {
		myimp = imp;
		if(!(arg==null||arg.equals(""))){
			File dest = new File(arg);
			directory = dest.getParent()+ System.getProperty("file.separator");
			name = dest.getName();
			System.out.println("saving as spider "+directory+name);
		}else{
			name=null;
		}
		return DOES_32 + NO_CHANGES;
	}


	/**
	 *  Main processing method for the SpiderWriter_ object
	 *
	 *@param  ip  Description of the Parameter
	 */
	public void run(ImageProcessor ip) {
		if(name==null){
			SaveDialog sd = new SaveDialog("save as...", myimp.getTitle(), ".spi");
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


		float[] head = createHeader();
		FileInfo fi = myimp.getFileInfo();
		fi.fileFormat = fi.RAW;
		//IJ.write("little endian=" + littleEndian);
		//fi.intelByteOrder = littleEndian;
		fi.fileName = name;
		fi.directory = directory;
		ImageWriter file = new ImageWriter(fi);
		try {
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(directory +System.getProperty("file.separator")+name)));
			for (float aHead : head) {
				out.writeFloat(aHead);
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
	private float[] createHeader() {
		int nsam = myimp.getWidth();
		int lenbyt = nsam * 4;
		int labrec = 1024 / lenbyt;
		if (1024 % lenbyt != 0) {
			labrec++;
		}
		int labbyt = labrec * lenbyt;

		float[] header = new float[labbyt / 4];
		header[0] = myimp.getNSlices();
		header[1] = myimp.getHeight();
		header[2] = header[1] + 1;
		header[3] = 1;
		header[4] = header[0] > 1 ? 3 : 1;
		header[11] = nsam;
		header[12] = labrec;
		header[21] = labbyt;
		header[22] = lenbyt;
		//header = computeMinMaxAvg();
		//setDateTime();
		//setTitle();
		// String tmp = "";
		// for (int i = 0; i < labbyt / 4; i++) {
		// tmp += head[i] + "\t";
		// if (i % 10 == 0) {
		// tmp += "\n";
		// }
		// }
		// System.out.println(tmp);
		return header;
	}

}
