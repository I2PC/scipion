/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package xmipp.ij;

import ij.IJ;

import java.awt.Menu;
import java.awt.MenuBar;
import java.awt.MenuItem;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.JCheckBoxMenuItem;

/**
 * 
 * @author Juanjo Vega
 */
public class XmippMenuBar extends MenuBar
{

	private Menu filemn;
	private Menu imagemn;
	private Menu advancedmn;
	private MenuItem savemi;
	private MenuItem saveasmi;
	private MenuItem openwithmi;
	private Menu infomn;
	private Menu adjustmn;
	private Menu transformmn;
	private Menu filtersmn;
	private Menu maskmn;
	private Menu thresholdingmn;
	private Menu binarymn;
	private Menu processmn;
	private Menu drawmn;
	private MenuItem propertiesmi;
	private MenuItem meanshiftmi;
	private MenuItem plotprofilemi;
	private MenuItem imagejmi;
	private MenuItem imageinfomi;
	private MenuItem fliphmi;
	private MenuItem flipvmi;
	private MenuItem cropmi;
	private MenuItem rotate90leftmi;
	private MenuItem rotate90rightmi;
	private ArrayList<String> requireij;
	private MenuItem bandpassmi;
	private MenuItem admi;

	public XmippMenuBar()
	{
		
		requireij = new ArrayList<String>();
		// menubar menus
		filemn = new Menu("File");
		imagemn = new Menu("Image");
		advancedmn = new Menu("Advanced");

		add(filemn);
		add(imagemn);
		add(advancedmn);

		// menubar file menu
		savemi = new MenuItem("Save");
		saveasmi = new MenuItem("Save As...");
		openwithmi = new MenuItem("Open with...");

		filemn.add(savemi);
		filemn.add(saveasmi);
		filemn.add(openwithmi);

		// menubar image menu
		infomn = new Menu("Info");
		adjustmn = new Menu("Adjust");
		transformmn = new Menu("Transform");
		filtersmn = new Menu("Filters");
		maskmn = new Menu("Mask");

		imagemn.add(infomn);
		imagemn.add(adjustmn);
		imagemn.add(transformmn);
		imagemn.add(filtersmn);
		imagemn.add(maskmn);

		// image info menu
		imageinfomi = new MenuItem("Image Info");
		addCommand(imageinfomi, "Show Info...");
		propertiesmi = new MenuItem("Properties");
		addCommand(propertiesmi, "Properties...");
		meanshiftmi = new MenuItem("Histogram");
		addCommand(meanshiftmi, "Histogram", true);// works only if imagej shown
		plotprofilemi = new MenuItem("Plot Profile");
		addCommand(plotprofilemi, "Plot Profile");// requires selection

		infomn.add(imageinfomi);
		infomn.add(propertiesmi);
		infomn.add(meanshiftmi);
		infomn.add(plotprofilemi);
		
		
		// image filters menu
		bandpassmi = new MenuItem("Bandpass Filter");
		addCommand(imageinfomi, "Bandpass Filter...");
		admi = new MenuItem("Anisotropic Diffussion");
		addCommand(admi, "Anisotropic Diffussion...");
		meanshiftmi = new MenuItem("Mean Shift");
		addCommand(meanshiftmi, "Mean Shift");

		filtersmn.add(bandpassmi);
		filtersmn.add(admi);
		filtersmn.add(meanshiftmi);
		
		

		// image transform menu
		cropmi = new MenuItem("Crop");
		addCommand(cropmi, "Crop", true);// works only if stack
		fliphmi = new MenuItem("Flip Horizontally");
		addCommand(fliphmi, "Flip Horizontally");
		flipvmi = new MenuItem("Flip Vertically");
		addCommand(flipvmi, "Flip Vertically");
		
		
		rotate90leftmi = new MenuItem("Rotate 90 Degrees Left");
		addCommand(rotate90leftmi, "Rotate 90 Degrees Left");
		rotate90rightmi = new MenuItem("Rotate 90 Degrees Right");
		addCommand(rotate90rightmi, "Rotate 90 Degrees Right");

		transformmn.add(cropmi);
		transformmn.add(fliphmi);
		transformmn.add(flipvmi);
		transformmn.add(rotate90leftmi);
		transformmn.add(rotate90rightmi);

		// menubar advanced menu
		imagejmi = new MenuItem("ImageJ");
		thresholdingmn = new Menu("Thresholding");
		binarymn = new Menu("Binary");
		processmn = new Menu("Process");
		drawmn = new Menu("Draw");

		advancedmn.add(imagejmi);
		advancedmn.add(thresholdingmn);
		advancedmn.add(binarymn);
		advancedmn.add(processmn);
		advancedmn.add(drawmn);
		
		
		//advanced threshold menu
		
		//advanced binary menu
		
		//advanced process menu
		
		//advanced drawn menu

		
		
		imagejmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				XmippImageWindow.openImageJ(Tool.VIEWER);
			}
		});

	}

	protected void addCommand(MenuItem mi, String command)
	{
		addCommand(mi, command, false);
	}
	
	
	protected void addCommand(MenuItem mi, String command, boolean isrequiredij)
	{
		if(isrequiredij)
			requireij.add(command);
		mi.setActionCommand(command);
		mi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				
				String command = ((MenuItem) e.getSource()).getActionCommand();
				if(requireij.contains(command))
					XmippImageWindow.openImageJ(Tool.VIEWER);
				if(command.equals("Anisotropic Diffusion..."))
					IJ.run("8-bit");
				IJ.run(command);
			}
		});
	}
	
	

}
