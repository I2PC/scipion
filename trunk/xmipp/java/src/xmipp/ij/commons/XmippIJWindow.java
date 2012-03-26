package xmipp.ij.commons;

import ij.ImagePlus;

public interface XmippIJWindow
{
	
	
	public void loadData();
	
	public void saveData() throws Exception;
	
	public void saveDataAs(String file) throws Exception;
	
	public ImagePlus getImagePlus();
	
	public String getImageFilePath();
	
}
