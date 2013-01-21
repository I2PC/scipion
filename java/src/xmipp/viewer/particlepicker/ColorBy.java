package xmipp.viewer.particlepicker;

import java.awt.Color;

import xmipp.jni.MetaData;

public class ColorBy
{
	int id;
	String name;
	Color color;
	private MetaData md;
	private double max;
	private double min;
	private double range;

	public ColorBy(int id, String name, Color color, MetaData md)
	{
		this.id = id;
		this.name = name;
		this.color = color;
		this.md = md;
		max = md.getColumnMax(id);
		min = md.getColumnMin(id);
		range = max - min;
	}
	
	public ColorBy(int id, String name, MetaData md)
	{
		this(id, name, getNextColor(), md);
	}


	private static Color[] colors = new Color[] { Color.BLUE, Color.CYAN, Color.GREEN, Color.MAGENTA, Color.ORANGE, Color.PINK, Color.YELLOW };
	private static int nextcolor;

	public static Color getNextColor()
	{
		Color next = colors[nextcolor];
		nextcolor++;
		if (nextcolor == colors.length)
			nextcolor = 0;
		return next;
	}
	
	public static Color[] getSampleColors() {
		return colors;
	}
	
	public String toString()
	{
		return name;
	}

	public Color getColor()
	{
		return color;
	}

	public void setColor(Color color2)
	{
		color = color2;
	}

	public String getName()
	{
		return name;
	}

	public int getId()
	{
		return id;
	}
	
	public Color getColor(double score)
	{
		
		double percent = 1;
		if(range > 0)
			percent = (score - min)/range;
		System.out.println(percent);
		Color scorecolor = new Color(color.getRed(), color.getGreen(), color.getBlue(), (int)(percent * 255));
		return scorecolor;
	}

}
