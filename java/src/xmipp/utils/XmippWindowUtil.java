/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *				Airen Zaldivar Peraza  (airenzp@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
package xmipp.utils;

import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JRootPane;

public class XmippWindowUtil
{

	/** Some colors contants */
	public static final Color LIGHT_BLUE = new Color(173, 216, 230);

	/**
	 * This function will be used to place the location of a windows relative to
	 * another windows or screen
	 */
	private static void setLocation(float positionx, float positiony, Container w, Dimension dim, Point offset)
	{
		Rectangle abounds = w.getBounds();
		int x = (int) (positionx * (dim.width - abounds.width) + offset.x);
		int y = (int) (positiony * (dim.height - abounds.height) + offset.y);
		w.setLocation(x, y);
	}

	public static void setLocation(float positionx, float positiony, Container w)
	{
		setLocation(positionx, positiony, w, w.getToolkit().getScreenSize(), new Point(0, 0));
	}

	public static void setLocation(float positionx, float positiony, Container w, Container parent)
	{
		setLocation(positionx, positiony, w, parent.getSize(), parent.getLocation());
	}

	/** Center the component in the screen */
	public static void centerWindows(Container w)
	{
		setLocation(0.5f, 0.5f, w);
	}

	/** Center the component relative to parent component */
	public static void centerWindows(Container w, Container parent)
	{
		setLocation(0.5f, 0.5f, w, parent.getSize(), parent.getLocation());
	}

	public static GridBagConstraints getConstraints(GridBagConstraints constraints, int x, int y)
	{
		return getConstraints(constraints, x, y, 1, 1);
	}

	public static GridBagConstraints getConstraints(GridBagConstraints constraints, int x, int y, int columns)
	{
		return getConstraints(constraints, x, y, columns, 1);
	}
	

	public static GridBagConstraints getConstraints(GridBagConstraints constraints, int x, int y, int columns, int rows, int fill)
	{
		constraints = getConstraints(constraints, x, y, columns, rows);
		constraints.fill = fill;
		return constraints;
	}
	

	public static JButton getIconButton(String icon, ActionListener listener){
		JButton btn = new JButton();
		btn.setIcon(XmippResource.getIcon(icon));
		Dimension dim = btn.getPreferredSize();
		dim.width = dim.height;
		btn.setPreferredSize(dim);
		btn.addActionListener(listener);
		return btn;
	}

	public static JButton getTextButton(String text, ActionListener listener)
	{
		JButton btn = new JButton(text);
		btn.setBackground(LIGHT_BLUE);
		btn.addActionListener(listener);
		return btn;
	}

	public static JLabel getIconLabel(String icon)
	{
		JLabel label = new JLabel();
		label.setIcon(XmippResource.getIcon(icon));
		return label;
	}

	public static GridBagConstraints getConstraints(GridBagConstraints constraints, int x, int y, int columns, int rows)
	{
		constraints.gridx = x;
		constraints.gridy = y;
		constraints.gridwidth = columns;
		constraints.gridheight = rows;
		if (columns > 1 && rows > 1)
		{
			constraints.fill = GridBagConstraints.BOTH;
		}
		else if (columns > 1)
		{
			constraints.fill = GridBagConstraints.HORIZONTAL;
		}
		else if (rows > 1)
		{
			constraints.fill = GridBagConstraints.VERTICAL;
		}
		else
		{
			constraints.fill = GridBagConstraints.NONE;
		}
		return constraints;
	}

	public static void openURI(String uri)
	{

		if (!java.awt.Desktop.isDesktopSupported())
			throw new IllegalArgumentException("Desktop is not supported (fatal)");

		if (uri == null)
			throw new IllegalArgumentException("Usage: OpenURI [URI [URI ... ]]");

		java.awt.Desktop desktop = java.awt.Desktop.getDesktop();

		if (!desktop.isSupported(java.awt.Desktop.Action.BROWSE))
		{

			throw new IllegalArgumentException("Desktop doesn't support the browse action (fatal)");
		}
		try
		{

			java.net.URI myuri = new java.net.URI(uri);
			desktop.browse(myuri);
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}
	}

	/** Block the gui and show the InfiniteProgressPanel */
	public static void blockGUI(JFrame window, String status)
	{
		final InfiniteProgressPanel progressPanel = new InfiniteProgressPanel(window, status);
		window.setGlassPane(progressPanel);
		progressPanel.start();
	}

	/** Release the gui from a previous block */
	public static void releaseGUI(JRootPane panel)
	{
		InfiniteProgressPanel progressPanel = (InfiniteProgressPanel) panel.getGlassPane();
		progressPanel.stop();
		progressPanel.setVisible(false);
	}

	

}