package xmipp.io;

import java.awt.*;
import java.io.*;
import ij.*;
import ij.io.*;
import ij.plugin.PlugIn;

public class SaveAsSPE_ implements PlugIn {

	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		if (imp == null) {
			IJ.showMessage("Save as SPE...", "No images are open.");
			return;
		}
		if (speType(imp.getFileInfo().fileType) < 0) {
			IJ.showMessage("Save as SPE...",
				"Supported types:\n" +
				"\n" +
				"32-bit Grayscale float : FLOAT\n" +
				"(32-bit Grayscale integer) : LONG\n" +
				"16-bit Grayscale integer: INT\n" +
				"(16-bit Grayscale unsigned integer) : UNINT\n");
			return;
		}
		String name = arg;
		if (arg == null || arg.equals("")) {
			name = imp.getTitle();
		}
		SaveDialog sd = new SaveDialog("Save as SPE...", name, ".spe");
		String file = sd.getFileName();
		if (file == null) return;
		String directory = sd.getDirectory();
		save(imp, directory, file);
	}

	public static void save(ImagePlus imp, String directory, String file) {
		if (imp == null) {
			IJ.showMessage("Save as SPE...", "No image selected.");
			return;
		}
		FileInfo fi = imp.getFileInfo();
		fi.intelByteOrder = true;
		int datatype = speType(fi.fileType);
		if (datatype < 0) {
			IJ.showMessage("Save as SPE...",
				"Supported types:\n" +
				"\n" +
				"32-bit Grayscale float : FLOAT\n" +
				"(32-bit Grayscale integer) : LONG\n" +
				"16-bit Grayscale integer: INT\n" +
				"(16-bit Grayscale unsigned integer) : UNINT\n");
			return;
		}
		SpeHeader header = new SpeHeader
			(datatype, imp.getWidth(), imp.getHeight(), imp.getStackSize());
		File f = new File(directory, file);
		try {
			FileOutputStream out = new FileOutputStream(f);
			byte[] h = header.getHeader();
			out.write(h, 0, h.length);
			ImageWriter writer = new ImageWriter(fi);
			writer.write(out);
			IJ.showStatus("");
		} catch (IOException e) {
			IJ.error("An error occured writing the file.\n \n" + e);
			IJ.showStatus("");
		}
	}

	public static int speType(int fiType) {
		switch (fiType) {
			case FileInfo.GRAY32_FLOAT:
				return SpeHeader.FLOAT;
			case FileInfo.GRAY32_INT:
				return SpeHeader.LONG;
			case FileInfo.GRAY16_SIGNED:
				return SpeHeader.INT;
			case FileInfo.GRAY16_UNSIGNED:
				return SpeHeader.UNINT;
			default:
				return -1;
		}
	}

}
