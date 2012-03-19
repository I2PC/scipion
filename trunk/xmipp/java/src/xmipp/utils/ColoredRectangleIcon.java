package xmipp.utils;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;

import javax.swing.Icon;


/**
 *  The ColoredRectangleIcon will paint a simple rectagle icon.
 *  The icon will ve painted with a specified color.
 *
 */
public class ColoredRectangleIcon implements Icon
{
	private Color color;
	private int width, height;
	private static final int gap = 3;

    public ColoredRectangleIcon(Color c, int w, int h)
    {
    	color = c;
    	width = w;
    	height = h;
    }//Constructor

    //  Implement the Icon Interface
	/**
	 *  Gets the width of this icon.
	 *
	 *  @return the width of the icon in pixels.
	 */
	@Override
    public int getIconWidth()
    {
		return width;
    }

	/**
	 *  Gets the height of this icon.
	 *
	 *  @return the height of the icon in pixels.
	 */
	@Override
    public int getIconHeight()
    {
		return height;
    }

   /**
    *  Paint the icons of this compound icon at the specified location
    *
    *  @param c The component on which the icon is painted
    *  @param g the graphics context
    *  @param x the X coordinate of the icon's top-left corner
    *  @param y the Y coordinate of the icon's top-left corner
    */
	@Override
    public void paintIcon(Component c, Graphics g, int x, int y)
    {
		x += gap;
		y += gap;
		int w = width - 2 * gap;
		int h = height - 2 * gap;
		g.setColor(color);
		g.fillRoundRect(x, y, w, h, 3, 3);
		g.setColor(Color.black);
		g.drawRoundRect(x, y, w, h, 3, 3);
		
    }

}//class ColoredRectangleIcon
