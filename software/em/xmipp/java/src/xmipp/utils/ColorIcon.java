package xmipp.utils;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;

import javax.swing.Icon;

/**
 * The ColoredRectangleIcon will paint a simple rectagle icon. The icon will ve
 * painted with a specified color.
 * 
 */
public class ColorIcon implements Icon
{
	private Color color;
	private int width, height;
	private int gap;
	private boolean rounded;
	private boolean withborder;
	private static final int ROUND_WIDTH = 3;

	public ColorIcon(Color c)
	{
		this(c, 16, 16);
	}

	public ColorIcon(Color c, int w, int h)
	{
		this(c, w, h, 0, true, true);
	}

	public ColorIcon(Color c, int w, int h, int g, boolean isRounded, boolean withborder)
	{
		color = c;
		width = w;
		height = h;
		gap = g;
		rounded = isRounded;
		this.withborder = withborder;

	}// Constructor

	// Implement the Icon Interface
	/**
	 * Gets the width of this icon.
	 * 
	 * @return the width of the icon in pixels.
	 */
	@Override
	public int getIconWidth()
	{
		return width;
	}

	/**
	 * Gets the height of this icon.
	 * 
	 * @return the height of the icon in pixels.
	 */
	@Override
	public int getIconHeight()
	{
		return height;
	}

	public Color getColor()
	{
		return color;
	}

	public void setColor(Color value)
	{
		color = value;
	}

	public int getGap()
	{
		return gap;
	}

	public void setGap(int g)
	{
		gap = g;
	}

	/**
	 * Paint the icons of this compound icon at the specified location
	 * 
	 * @param c
	 *            The component on which the icon is painted
	 * @param g
	 *            the graphics context
	 * @param x
	 *            the X coordinate of the icon's top-left corner
	 * @param y
	 *            the Y coordinate of the icon's top-left corner
	 */
	@Override
	public void paintIcon(Component c, Graphics g, int x, int y)
	{
		x += gap;
		y += gap;
		int w = width - 2 * gap;
		int h = height - 2 * gap;
		g.setColor(color);
		if (rounded)
			g.fillRoundRect(x, y, w, h, ROUND_WIDTH, ROUND_WIDTH);
		else
			g.fillRect(x, y, w, h);
		if (withborder)
		{
			g.setColor(Color.black);
			if (rounded)
				g.drawRoundRect(x, y, w, h, ROUND_WIDTH, ROUND_WIDTH);
			else
				g.drawRect(x, y, w, h);
		}
	}

}// class ColoredRectangleIcon
