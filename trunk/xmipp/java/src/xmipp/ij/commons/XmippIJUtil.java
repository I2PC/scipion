package xmipp.ij.commons;

import ij.IJ;
import xmipp.jni.Filename;

public class XmippIJUtil {
	
	
	private static XmippImageJ xij;
	
	public static XmippImageJ showImageJ(Tool tool)
	{
		if (IJ.getInstance() == null)
		{
			xij = new XmippImageJ();
			IJ.run("Install...", "install=" + Filename.getXmippPath("java/src/xmipp/ij/commons/XmippMacros.txt"));
		}
		else if (!xij.isVisible())
			xij.setVisible(true);
		return xij;
	}
	
	public static XmippImageJ getXmippImageJ()
	{
		return xij;
	}


}
