package xmipp.ij.commons;

public enum Tool
{
	IMAGEJ, PICKER, VIEWER;
	
	
	public static Tool getTool(String name)
	{
		if (name.equalsIgnoreCase("Particle Picker Tool"))
			return Tool.PICKER;
		else if(name.equalsIgnoreCase("Xmipp Micrograph Viewer Tool"))
			return Tool.VIEWER;
		return Tool.IMAGEJ;
	}
	
	public static String getTool(Tool tool)
	{
		if(tool == Tool.PICKER)
			return "Particle Picker Tool";
		else if(tool == Tool.VIEWER)
			return "Xmipp Micrograph Viewer Tool";
		return null;
	}
}
