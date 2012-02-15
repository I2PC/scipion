package xmipp.utils;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.Point;
import java.awt.Rectangle;


public class WindowUtils {

	/** This function will be used to place the location of 
	 * a windows relative to another windows or screen
	 */
	private static void setLocation(float positionx, float positiony, Container w, Dimension dim, Point offset) {
		Rectangle abounds = w.getBounds();
		int x = (int) (positionx * (dim.width - abounds.width) + offset.x);
		int y = (int) (positiony * (dim.height - abounds.height) + offset.y);
		w.setLocation(x, y);
	}

	public static void setLocation(float positionx, float positiony, Container w) {
		setLocation(positionx, positiony, w, w.getToolkit().getScreenSize(), new Point(0, 0));
	}
	
	/** Center the component in the screen */
	public static void centerWindows(Container w){
		setLocation(0.5f, 0.5f, w);		
	}
	
	/** Center the component relative to parent component */
	public static void centerWindows(Container w, Container parent){
		setLocation(0.5f, 0.5f, w, parent.getSize(), parent.getLocation());		
	}	

	public static GridBagConstraints getConstraints(
			GridBagConstraints constraints, int x, int y, int columns) {
		constraints.gridx = x;
		constraints.gridy = y;
		constraints.gridwidth = columns;
		constraints.fill = (columns > 1) ? GridBagConstraints.HORIZONTAL
				: GridBagConstraints.NONE;
		return constraints;
	}

	public static void openURI(String uri) {

		if (!java.awt.Desktop.isDesktopSupported())
			throw new IllegalArgumentException(
					"Desktop is not supported (fatal)");

		if (uri == null)
			throw new IllegalArgumentException(
					"Usage: OpenURI [URI [URI ... ]]");

		java.awt.Desktop desktop = java.awt.Desktop.getDesktop();

		if (!desktop.isSupported(java.awt.Desktop.Action.BROWSE)) {

			throw new IllegalArgumentException("Desktop doesn't support the browse action (fatal)");
		}
		try {

			java.net.URI myuri = new java.net.URI(uri);
			desktop.browse(myuri);
		} catch (Exception e) {
			throw new IllegalArgumentException(e);
		}
	}
}
