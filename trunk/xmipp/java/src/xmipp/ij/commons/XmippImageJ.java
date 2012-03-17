package xmipp.ij.commons;

import ij.Executer;

import java.awt.event.WindowEvent;

import ij.ImageJ;

public class XmippImageJ extends ImageJ
{
	public void windowClosing(WindowEvent e) {
		setVisible(false);
	}
	
	public void close()
	{
		new Executer("Quit", null);
		
	}
}
