import java.awt.*;
import java.awt.event.*;
import java.io.*;
import ij.*;
import ij.io.*;
import ij.plugin.PlugIn;

public class OpenSPE_ implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open SPE...", arg);
		String file = od.getFileName();
		if (file == null) return;
		String directory = od.getDirectory();
		ImagePlus imp = open(directory, file);
		if (imp != null ) {
			imp.show();
		} else {
			IJ.showMessage("Open SPE...", "Failed.");
		}
	}

	public static ImagePlus open(String directory, String file) {
		File f = new File(directory, file);
		try {
			FileInputStream in = new FileInputStream(f);
			byte[] h = new byte[SpeHeader.headerSize];
			in.read(h, 0, h.length);
			SpeHeader header = new SpeHeader(h);
			int speType = header.getDatatype();
			int fiType = fileInfoType(speType);
			if (fiType < 0) {
				IJ.showMessage("Open SPE...",
					"Invalid data type.");
				return null;
			}

			//if (speType == header.UNINT) {
			//	boolean convert = IJ.showMessageWithCancel("Open SPE...",
			//		"Convert UNSIGNED 16-bit integer \n to SIGNED 16-bit integer ?");
			//	if(convert)
			//	fiType = FileInfo.GRAY16_SIGNED;
			//}

			FileInfo fi = new FileInfo();
			fi.directory = directory;
			fi.fileFormat = fi.RAW;
			fi.fileName = file;
			fi.fileType = fiType;
			fi.gapBetweenImages = 0;
			fi.height = header.getHeight();
			fi.intelByteOrder = true;
			fi.nImages = header.getStackSize();
			fi.offset = SpeHeader.RAW_OFFSET;
			fi.width = header.getWidth();
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
