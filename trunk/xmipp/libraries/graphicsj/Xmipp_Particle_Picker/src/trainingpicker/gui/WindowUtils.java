package trainingpicker.gui;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.Rectangle;
import java.awt.Window;
import java.net.URI;
import java.util.logging.Level;
import java.awt.Desktop;

import trainingpicker.model.TrainingPicker;


public class WindowUtils {

	public static void centerScreen(double positionx, double positiony, Container w) {
		Dimension dim = w.getToolkit().getScreenSize();
		Rectangle abounds = w.getBounds();
		int x = (int) (positionx * (dim.width - abounds.width));
		int y = (int) (positiony * (dim.height - abounds.height));
		w.setLocation(x, y);

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
