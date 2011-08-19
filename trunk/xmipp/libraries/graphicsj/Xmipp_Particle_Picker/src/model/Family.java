package model;
import java.awt.Color;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;


public class Family {
	
	private String name;
	private Color color;
	private int size;
	
	private static int sizemax = 1000;
	private static Family dfamily = new Family("Any", Color.green);
	private static Color[] colors = new Color[]{Color.BLUE, Color.CYAN, 
										Color.GREEN,
										Color.MAGENTA, Color.ORANGE, 
										Color.PINK, Color.RED, Color.YELLOW};
	
	
	public Family(String name, Color color, int size)
	{
		if(size > sizemax)
			throw new IllegalArgumentException(String.format("Max size is %s, %s not allowed", sizemax, size));
		if (name == null || name.equals(""))
			throw new IllegalArgumentException(Constants.getEmptyFieldMsg("name"));
		this.name = name;
		this.color = color;
		this.size = size;
	}
	
	public Family(String name, Color color)
	{
		this.name = name;
		this.color = color;
		this.size = getDefaultSize();
		
	}
	
	public static String getOFilename()
	{
		return PPConfiguration.getOutputPath("families.xmd");
	}
	
	public int getSize() {
		return size;
	}


	public void setSize(int size) {
		if(size > sizemax)
			throw new IllegalArgumentException(String.format("Max size is %s, %s not allowed", sizemax, size));
		this.size = size;
	}


	public static Family getDefaultgp() {
		return dfamily;
	}
	

	public String getName() {
		return name;
	}

	public void setName(String name) {
		if (name == null || name.equals(""))
			throw new IllegalArgumentException(Constants.getEmptyFieldMsg("name"));
		this.name = name;
	}

	public Color getColor() {
		return color;
	}

	public void setColor(Color color) {
		this.color = color;
	}
	
	public static Family getDefaultFamily()
	{
		return dfamily;
	}
	
	public String toString()
	{
		return name;
	}
	
	public static Color[] getSampleColors()
	{
		return colors;
	}
	
	public static Color getColor(String name)
	{
		Color color;
		try {
		    Field field = Class.forName("java.awt.Color").getField(name);
		    color = (Color)field.get(null);//static field, null for parameter
		} catch (Exception e) {
		    color = null; // Not defined
		}
		return color;
	}
	
	public static int getDefaultSize()
	{
		return 100;
	}
}
