package xmipp.ij.commons;


public interface XmippIJWindow
{
	public void loadData();
	
	public void saveData() throws Exception;
	
	public void saveDataAs(String file) throws Exception;
	
	public ImagePlusLoader getImagePlusLoader();
	
	public boolean isVolume();

	public boolean isStack();
	
	public void openMaskToolbar();

	public XmippImageCanvas getCanvas();
	
}
