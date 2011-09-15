package trainingpicker.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;

import javax.swing.Icon;

public class ColorIcon implements Icon {
	
	private static int width = 14;
	private static int height = 14;
	private Color color;
	
	public ColorIcon(Color color)
	{
		this.color = color;
	}

	@Override
	public int getIconHeight() {
		return height;
	}

	@Override
	public int getIconWidth() {
		return width;
	}

	@Override
	public void paintIcon(Component c, Graphics g, int x, int y) {
		 g.setColor(color);  
	     g.fillRect(x, y, width, height);  
	     g.setColor(Color.black);  
	     g.drawRect(x, y, width, height); 

	}
	
	

}
