package gui;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.Rectangle;
import java.awt.Window;

public class WindowUtils {
	
	public static void centerScreen(double position, Container w) {
		Dimension dim = w.getToolkit().getScreenSize();
		Rectangle abounds = w.getBounds();
		int x = (int)(position * (dim.width - abounds.width)) ;
		int y = (dim.height - abounds.height) / 2;
		w.setLocation(x, y);
		
	}
	
	public static GridBagConstraints updateConstraints(GridBagConstraints constraints, int x, int y, int columns)
	{
		constraints.gridx = x;
		constraints.gridy = y;
		constraints.gridwidth = columns;
		constraints.fill = (columns > 1)? GridBagConstraints.HORIZONTAL: GridBagConstraints.NONE;
		return constraints;
	}

}
