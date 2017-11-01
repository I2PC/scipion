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
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.awt.GridBagConstraints;
import java.awt.Image;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Window;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.net.UnknownHostException;
import java.util.Arrays;

import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JRootPane;

public class XmippWindowUtil
{
        public static final Color firebrick = Color.decode("#B22222");
        public static final Color lightgrey = Color.decode("#EAEBEC");

	/** Some colors contants */
	public static final Color LIGHT_BLUE = new Color(173, 216, 230);
        
        

	/**
	 * This function will be used to place the location of a windows relative to
	 * another windows or screen
	 */
	private static void setLocation(double positionx, double positiony, Container w, Rectangle bounds)
	{
		Rectangle abounds = w.getBounds();
		int x = (int) (positionx * (bounds.width - abounds.width) + bounds.x);
		int y = (int) (positiony * (bounds.height - abounds.height) + bounds.y);
		w.setLocation(x, y);
	}

	public static void setLocation(double positionx, double positiony, Window w)
	{
        setLocation(positionx, positiony, w, getScreenRectangle());
	}
	
	public static Rectangle getScreenRectangle()
	{
		GraphicsDevice gd = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice();
		return gd.getDefaultConfiguration().getBounds();
	}

	public static void setLocation(float positionx, float positiony, Container w, Container parent)
	{
		setLocation(positionx, positiony, w, parent.getBounds());
	}

	/** Center the component in the screen */
	public static void centerWindows(Window w)
	{
		setLocation(0.5f, 0.5f, w);
	}

	/** Center the component relative to parent component */
	public static void centerWindows(Container w, Container parent)
	{
		setLocation(0.5f, 0.5f, w, parent.getBounds());
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
	
	public static JButton getTextIconButton(String text, String icon, ActionListener listener){
		JButton btn = new JButton(text, XmippResource.getIcon(icon));
		btn.addActionListener(listener);
		return btn;
	}

	public static JButton getTextButton(String text, ActionListener listener)
	{
		JButton btn = new JButton(text);
                //if(!isScipion)
                //    btn.setBackground(LIGHT_BLUE);
		btn.addActionListener(listener);
		return btn;
	}

	public static JLabel getIconLabel(String icon)
	{
		JLabel label = new JLabel();
		label.setIcon(XmippResource.getIcon(icon));
		return label;
	}
        
    public static JButton getScipionIconButton(String text) {
        Icon icon = XmippResource.getIcon("fa-plus-circle.png");
        JButton button = new JButton(text.replace("Create ", ""), icon);
        button.setToolTipText(text);
        button.setBackground(firebrick);
        button.setForeground(Color.WHITE);
        return button;
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

    /** Block the gui and show the InfiniteProgressPanel without hiding the content */
    public static void blockGUI(JFrame window, String status){

        blockGUI(window, status, false);
    }

    /** Block the gui and show the InfiniteProgressPanel */
    /**
     * The hide content parameter is needed for the galleryData export process.
     *  If the content is hidden, there will not be any refresh and therefore there will not be a concurrency access to
     *  the metadata object in different threads (painting the window and exporting).
     *
    **/
	public static void blockGUI(JFrame window, String status, boolean hideContent)
	{
		final InfiniteProgressPanel progressPanel = new InfiniteProgressPanel(window, status);
		window.setGlassPane(progressPanel);
        if (hideContent) {
            window.getLayeredPane().setVisible(false);
        }
		progressPanel.start();
	}

	/** Release the gui from a previous block */
	public static void releaseGUI(JRootPane panel)
	{
		InfiniteProgressPanel progressPanel = (InfiniteProgressPanel) panel.getGlassPane();
		progressPanel.stop();
		progressPanel.setVisible(false);
        panel.getLayeredPane().setVisible(true);
	}

    public static void executeGUICommand(final String[] command, final JFrame frame, String msg)
    {


        XmippWindowUtil.blockGUI(frame, msg);
        new Thread(new Runnable() {

            @Override
            public void run() {

                try {

                    String output = executeCommand(command, true);
                    XmippWindowUtil.releaseGUI(frame.getRootPane());

                    if(output != null && !output.isEmpty())
                        System.out.println(output);

                } catch (Exception ex) {
                    throw new IllegalArgumentException(ex.getMessage());
                }

            }
        }).start();
    }
        
    public static String executeCommand(String[] command, boolean wait) throws Exception {
            //System.out.println(Arrays.toString(command));
            Process p = Runtime.getRuntime().exec(command);
            if(wait)
            {
                p.waitFor();
                return readProcessOutput(p);
            }
            return null;
    }
    
    public static String executeCommand(String command, boolean wait, String dir) throws Exception {
        //System.out.println(command);
        Process p = Runtime.getRuntime().exec(command, null, new File(dir));
        if(wait)
        {
            p.waitFor();
            return readProcessOutput(p);
        }
        return null;
}
    
    public static String executeCommand(String command, boolean wait) throws Exception {
            //System.out.println(command);
            Process p = Runtime.getRuntime().exec(new String[] { "bash", "-c", command });
            if(wait)
            {
                p.waitFor();
                return readProcessOutput(p);
            }
            return null;
    }
        
    public static void runCommand(String command, Integer port) {
        //System.out.println(command);
    	if(port == null)
    	{
    		XmippDialog.showError(null, XmippMessage.getEmptyFieldMsg("Scipion port"));
    		return;
    	}
        String hostName = "";
 
        try {
            Socket socket = new Socket(hostName, port);
            PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
            BufferedReader in = new BufferedReader(new InputStreamReader(socket.getInputStream()));
            out.println(command);
            in.readLine();
            socket.close();
        
            
        } catch (UnknownHostException e) {
            System.err.println("Don't know about host " + hostName);
        } catch (IOException e) {
            System.err.println("Couldn't get I/O for the connection to " +
                hostName);
        }
        }



    public static String readProcessOutput(Process p) throws IOException
    {
        StringBuffer output = new StringBuffer();
        BufferedReader reader
                = new BufferedReader(new InputStreamReader(p.getInputStream()));


        String line = "";
        while ((line = reader.readLine()) != null) {
            output.append(line + "\n");
        }
        reader = new BufferedReader(new InputStreamReader(p.getErrorStream()));

        while ((line = reader.readLine()) != null) {
            output.append(line + "\n");
        }
        return output.toString();

    }
    
    public static void setScipionImageIcon(Window w)
    {
            Image img = XmippResource.getIcon("scipion_logo.png").getImage();
            w.setIconImage(img);
    }

}