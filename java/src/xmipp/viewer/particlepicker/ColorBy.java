package xmipp.viewer.particlepicker;

import java.awt.Color;

public class ColorBy
{
	int id;
	String name;
	Color color;

	public ColorBy(int id, String name, Color color)
	{
		this.id = id;
		this.name = name;
		this.color = color;
	}
	
	public ColorBy(int id, String name)
	{
		this(id, name, getNextColor());
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
	
	public static Color getColor(double score, Color color)
	{
		Color scorecolor = new Color(color.getRed(), color.getGreen(), color.getBlue(), (int)(Math.min(1, score) * 255));
		return scorecolor;
	}

}
