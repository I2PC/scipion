package xmipp.ij.commons;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.ImageWindow;
import ij.io.FileSaver;
import java.io.File;
import java.util.Date;
import xmipp.jni.Filename;

@SuppressWarnings("serial")
public class XmippImageWindow extends ImageWindow implements XmippIJWindow {

	public static void main(String[] args) {
		try {
			// openImageJ(Tool.VIEWER);
			XmippImageWindow w = new XmippImageWindow(
					"/home/airen/Coss/Xmipp/BPV_2/InputData/BPV_1386.mrc");
//			IJ.open( "/home/airen/Coss/Xmipp/BPV_2/InputData/BPV_1386.mrc");

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private String file;
	private long modified;

	public XmippImageWindow(String file) throws Exception {
		this(XmippImageConverter.loadImage(file), file);
		this.file = file;
		this.modified = new File(file).lastModified();
	}

	public XmippImageWindow(ImagePlus imp) {
		this(imp, "");
	}

	public XmippImageWindow(ImagePlus imp, String title) {
		super(imp, new XmippImageCanvas(imp));
		// file = imp.getOriginalFileInfo().directory + File.separator +
		// imp.getOriginalFileInfo().fileName;
		setTitle(title);
		setMenuBar(new XmippMenuBar(this));
	}

	public static void openImageJ(Tool tool) {
		if (IJ.getInstance() == null) {
			new ImageJ();
			IJ.run("Install...",
					"install="
							+ Filename
									.getXmippPath("java/src/xmipp/ij/commons/XmippMacros.txt"));
		}
		boolean recognized = IJ.setTool(Tool.getTool(tool));
		// System.out.println(recognized);
	}

	@Override
	public void loadData() {
		try {
			if (file != null && new File(file).lastModified() > modified) {
				ImagePlus imp = XmippImageConverter.loadImage(file);
				setImage(imp);// second alone does not work
				updateImage(imp);// first one alone does not work
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	@Override
	public void saveDataAs(String file) throws Exception {
		this.file = file;
		XmippImageConverter.writeImagePlus(imp, file);
	}

	@Override
	public void saveData() throws Exception {
		saveDataAs(imp.getTitle());
	}
}// class XmippImageWindow
