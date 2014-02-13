/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package xmipp.ij.commons;

import ij.IJ;
import ij.ImagePlus;
import ij.process.StackConverter;
import ij3d.Content;
import ij3d.Image3DUniverse;

import java.awt.CheckboxMenuItem;
import java.awt.Frame;
import java.awt.Menu;
import java.awt.MenuBar;
import java.awt.MenuItem;
import java.awt.MenuShortcut;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.io.File;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;

import javax.swing.JOptionPane;
import javax.vecmath.Color3f;

import xmipp.utils.QuickHelpJDialog;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;

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
	private MenuItem openwith3dmi;
	private Menu infomn;
	private Menu adjustmn;
	private Menu transformmn;
	private Menu filtersmn;
	private Menu thresholdingmn;
	private Menu binarymn;
	private Menu processmn;
	private Menu drawmn;
	private MenuItem imagejmi;
	private Menu profilemn;
	private CheckboxMenuItem pollmi;
	protected Object timer;
	private final XmippIJWindow xw;
	private PollTimer polltimer;
	private MenuItem refreshmi;
	private CheckboxMenuItem wrapmi;
	private CheckboxMenuItem ugmi;
	private MenuItem exitmi;
	private Menu helpmn;
	private MenuItem keyassistmi;
	private QuickHelpJDialog keyassistdlg;

	enum IJRequirement
	{
		BINARY, EIGHTBIT, IMAGEJ, STACK, THIRTYTWOBIT, RGB, VOLUME

	};

	public XmippMenuBar(XmippIJWindow xw)
	{
		this.xw = xw;
		// menubar menus
		filemn = new Menu("File");
		imagemn = new Menu("Image");
		advancedmn = new Menu("Advanced");
		helpmn = new Menu("Help");

		add(filemn);
		add(imagemn);
		add(advancedmn);
		add(helpmn);

		// menubar file menu
		savemi = new MenuItem("Save");
		savemi.setShortcut(new MenuShortcut(KeyEvent.VK_S));
		savemi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				try
				{
					if(XmippMenuBar.this.xw.getImagePlusLoader().existsFile())
						XmippMenuBar.this.xw.saveData();
					else
						XmippMenuBar.this.saveAs();
				}
				catch (Exception ex)
				{
					ex.printStackTrace();
					XmippDialog.showInfo(null, ex.getMessage());
				}

			}
		});
		saveasmi = new MenuItem("Save As...");
		saveasmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				saveAs();

			}
		});

		openwith3dmi = new MenuItem("Open with 3D Viewer");
		openwith3dmi.setEnabled(xw.isVolume());
//		openwith3dmi.setEnabled(Filename.isVolume(xw.getImageFilePath()));
		openwith3dmi.addActionListener(new ActionListener()
		{
			
			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				ImagePlus imp = XmippMenuBar.this.xw.getImagePlusLoader().getImagePlus();
				if (imp.getImageStackSize() == 1)
					XmippDialog.showInfo(null, "Only for Stack");
				else
					openImagePlusAs3D(imp);
				
			}
		});
		filemn.add(openwith3dmi);
		addIJMenuItem(filemn, "Open with Volume Viewer/3D Slicer", "Volume Viewer", IJRequirement.VOLUME);
		addIJMenuItem(filemn, "Open with VolumeJ", "VolumeJ ", IJRequirement.VOLUME);
		refreshmi = new MenuItem("Refresh");
		refreshmi.setEnabled(xw.getImagePlusLoader().allowsPoll());
		refreshmi.setShortcut(new MenuShortcut(KeyEvent.VK_F5));
		refreshmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				XmippMenuBar.this.xw.loadData();

			}
		});

		pollmi = new CheckboxMenuItem("Poll");
		pollmi.setEnabled(xw.getImagePlusLoader().allowsPoll());
		pollmi.addItemListener(new ItemListener()
		{

			@Override
			public void itemStateChanged(ItemEvent e)
			{
				boolean poll = pollmi.getState();
				if (poll)
					poll();
				else
					polltimer.stop();
			}
		});
		
		ugmi = new CheckboxMenuItem("Use Geometry");
		ugmi.setEnabled(xw.getImagePlusLoader().allowsGeometry());
		ugmi.setState(xw.getImagePlusLoader().getUseGeometry());
		ugmi.addItemListener(new ItemListener()
		{

			@Override
			public void itemStateChanged(ItemEvent e)
			{
				boolean ug = ugmi.getState();
				useGeometry(ug);
				
			}
		});

		
		wrapmi = new CheckboxMenuItem("Wrap");
		wrapmi.setEnabled(xw.getImagePlusLoader().allowsGeometry());
		wrapmi.setState(xw.getImagePlusLoader().isWrap());
		wrapmi.addItemListener(new ItemListener()
		{

			@Override
			public void itemStateChanged(ItemEvent e)
			{
				wrap(wrapmi.getState());
				
			}
		});
		
		exitmi = new MenuItem("Exit");
		exitmi.setShortcut(new MenuShortcut(KeyEvent.VK_Q));
		exitmi.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {

				close();
			}
		});



		filemn.add(savemi);
		filemn.add(saveasmi);
		addIJMenuItem(filemn, "Duplicate", "Duplicate...", IJRequirement.IMAGEJ);
		filemn.add(refreshmi);
		filemn.add(pollmi);
		filemn.add(ugmi);
		filemn.add(wrapmi);
		filemn.add(exitmi);

		// menubar image menu
		infomn = new Menu("Info");
		adjustmn = new Menu("Adjust");
		transformmn = new Menu("Transform");
		filtersmn = new Menu("Filters");

		imagemn.add(infomn);
		imagemn.add(adjustmn);
		imagemn.add(transformmn);
		imagemn.add(filtersmn);
		addIJMenuItem(imagemn, "Masks Tool Bar", "Masks Tool Bar", IJRequirement.IMAGEJ);// missing
																							// plugin

		// image info menu
		addIJMenuItem(infomn, "Show Info", "Show Info...");
		addIJMenuItem(infomn, "Properties", "Properties...");
		addIJMenuItem(infomn, "Histogram", "Histogram", IJRequirement.IMAGEJ);
		addIJMenuItem(infomn, "Plot Profile", "Plot Profile", IJRequirement.IMAGEJ);

		// image adjust menu
		addIJMenuItem(adjustmn, "Brightness/Contrast", "Brightness/Contrast...");
		addIJMenuItem(adjustmn, "Enhance Contrast", "Enhance Contrast");

		addIJMenuItem(adjustmn, "Crop", "Crop", IJRequirement.IMAGEJ);
		addIJMenuItem(adjustmn, "Scale", "Scale...");
		addIJMenuItem(adjustmn, "Untilt Stack", "Untilt Stack", IJRequirement.STACK);

		addIJMenuItem(adjustmn, "Reslice", "Reslice [/]...", IJRequirement.VOLUME);

		// image transform menu
		addIJMenuItem(transformmn, "Flip Horizontally", "Flip Horizontally");
		addIJMenuItem(transformmn, "Flip Vertically", "Flip Vertically");
		addIJMenuItem(transformmn, "Rotate 90 Degrees Left", "Rotate 90 Degrees Left");
		addIJMenuItem(transformmn, "Rotate 90 Degrees Right", "Rotate 90 Degrees Right");

		// image filters menu
		addIJMenuItem(filtersmn, "Bandpass Filter", "Bandpass Filter...");
		addIJMenuItem(filtersmn, "Anisotropic Diffusion", "Anisotropic Diffusion...", IJRequirement.EIGHTBIT);
		addIJMenuItem(filtersmn, "Mean Shift", "Mean Shift");

		// menubar advanced menu
		imagejmi = new MenuItem("ImageJ");
		imagejmi.setShortcut(new MenuShortcut(KeyEvent.VK_I));
		thresholdingmn = new Menu("Thresholding");
		binarymn = new Menu("Binary");
		processmn = new Menu("Process");
		drawmn = new Menu("Draw");
		profilemn = new Menu("Profile");

		advancedmn.add(imagejmi);
		advancedmn.add(thresholdingmn);
		advancedmn.add(binarymn);
		advancedmn.add(processmn);
		advancedmn.add(drawmn);
		advancedmn.add(profilemn);
		
		
		keyassistmi = new MenuItem("Key Assist");
		keyassistmi.addActionListener(new ActionListener()
		{
			
			

			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				try
				{
				if(keyassistdlg == null)
					keyassistdlg = new QuickHelpJDialog(null, false, "Key Assist...", getKeyAssist());
				keyassistdlg.setVisible(true);
				}
				catch(Exception e)
				{
					XmippDialog.showInfo(null, e.getMessage());
				}
			}
		});
		helpmn.add(keyassistmi);
		

		// advanced threshold menu
		addIJMenuItem(thresholdingmn, "Threshold", "Threshold...");
		addIJMenuItem(thresholdingmn, "Otsu Threshold", "Otsu Thresholding", IJRequirement.EIGHTBIT);
		addIJMenuItem(thresholdingmn, "Multi Otsu Threshold", "Multi OtsuThreshold");
		addIJMenuItem(thresholdingmn, "Maximum Entropy Threshold", "Entropy Threshold", IJRequirement.EIGHTBIT);
		addIJMenuItem(thresholdingmn, "Mixture Modeling Threshold", "Mixture Modeling", IJRequirement.EIGHTBIT);
		addIJMenuItem(thresholdingmn, "Robust Automatic Threshold Selection", "RATS ");
		//addIJMenuItem(thresholdingmn, "Simple Iterative Object Extraction", "SIOX Segmentation", IJRequirement.RGB);

		// advanced binary menu
		addIJMenuItem(binarymn, "Voxel Counter", "Voxel Counter", IJRequirement.STACK);// stack
																						// required
		addIJMenuItem(binarymn, "Erode", "Erode", IJRequirement.BINARY);
		addIJMenuItem(binarymn, "Dilate", "Dilate", IJRequirement.BINARY);
		addIJMenuItem(binarymn, "Open", "Open", IJRequirement.BINARY);
		addIJMenuItem(binarymn, "Close", "Close", IJRequirement.BINARY);
		addIJMenuItem(binarymn, "Float Morphology", "Float Morphology", IJRequirement.THIRTYTWOBIT);// missing
																									// plugin
		addIJMenuItem(binarymn, "Outline", "Outline", IJRequirement.BINARY, IJRequirement.EIGHTBIT);
		addIJMenuItem(binarymn, "Fill Holes", "Fill Holes", IJRequirement.BINARY, IJRequirement.EIGHTBIT);
		addIJMenuItem(binarymn, "Skeletonize", "Skeletonize", IJRequirement.BINARY, IJRequirement.EIGHTBIT);
		addIJMenuItem(binarymn, "Distance Map", "Distance Map", IJRequirement.BINARY, IJRequirement.EIGHTBIT);
		addIJMenuItem(binarymn, "Ultimate Points", "Ultimate Points", IJRequirement.BINARY, IJRequirement.EIGHTBIT);
		addIJMenuItem(binarymn, "Watershed", "Watershed");
		addIJMenuItem(binarymn, "Voronoi", "Voronoi", IJRequirement.BINARY, IJRequirement.EIGHTBIT);

		// advanced process menu

		addIJMenuItem(processmn, "Subtract Background", "Subtract Background...");
		addIJMenuItem(processmn, "Gaussian Blur", "Gaussian Blur...");
		addIJMenuItem(processmn, "Convolve", "Convolve...");
		addIJMenuItem(processmn, "Median", "Median...");
		addIJMenuItem(processmn, "FFT", "FFT");// memory error
		addIJMenuItem(processmn, "Straighten Curved Objects", "Straighten...", IJRequirement.IMAGEJ);

		// advanced drawn menu
		addIJMenuItem(drawmn, "Dotted Line", "Dotted Line", IJRequirement.IMAGEJ);
		addIJMenuItem(drawmn, "Radial Grid", "Radial Grid");
		addIJMenuItem(drawmn, "Concenctric Circles", "Concentric Circles");

		// advanced profile menu
		addIJMenuItem(profilemn, "Line Analyzer", "Line Analyzer", IJRequirement.IMAGEJ);
		addIJMenuItem(profilemn, "Oval Profile Plot", "Oval Profile", IJRequirement.IMAGEJ);
		addIJMenuItem(profilemn, "Radial Profile Plot Angle", "Radial Profile Angle", IJRequirement.IMAGEJ);
		addIJMenuItem(profilemn, "Radial Profile Plot Height", "Radial Profile Height", IJRequirement.IMAGEJ);
		addIJMenuItem(profilemn, "Contour Plotter", "ContourPlotter ");

		imagejmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				XmippUtil.showImageJ(Tool.VIEWER);
			}
		});

	}
        

        protected void close()
        {
            Frame w = (Frame)xw;
            w.setVisible(false);
            w.dispose();
            XmippApplication.removeInstance();
        }

	
	protected void useGeometry(boolean ug)
	{
		xw.getImagePlusLoader().setUseGeometry(ug);
		xw.loadData();
	}

	private void saveAs()
	{
		XmippFileChooser fc = new XmippFileChooser();
		int returnVal = fc.showOpenDialog(null);

		try
		{
			if (returnVal == XmippFileChooser.APPROVE_OPTION)
			{
				File file = fc.getSelectedFile();
				XmippMenuBar.this.xw.saveDataAs(file.getAbsolutePath());
			}
		}
		catch (Exception ex)
		{
			XmippDialog.showInfo(null, ex.getMessage());
		}
	}
	


	protected void wrap(boolean value)
	{
		xw.getImagePlusLoader().setWrap(value);
		xw.loadData();
	}
	
	public static void openImagePlusAs3D(ImagePlus ip) {
		try {
			int UNIVERSE_W = 400, UNIVERSE_H = 400;
			// Checks if java3D is available or not.
			Class.forName("javax.media.j3d.J3DBuffer");

			new StackConverter(ip).convertToRGB();

			Image3DUniverse universe = new Image3DUniverse(UNIVERSE_W,
					UNIVERSE_H);

			// Adds the sphere image plus to universe.
			Content c = universe.addSurfacePlot(ip, new Color3f(1f, 165f / 255,
					82f / 255), "1", 50, new boolean[] { true, true, true }, 1);
			c.displayAs(Content.SURFACE);
			c.setColor(new Color3f(1f, 165f / 255, 82f / 255));

			universe.show(); // Shows...
		} catch (final ClassNotFoundException e) {
			IJ.error("Java 3D not found. Please, check your installation.");
		}
	}

	private void addIJMenuItem(Menu mn, String name, String command, IJRequirement... requirements)
	{
		MenuItem mi = new MenuItem(name);
		addCommand(mi, command, requirements);
		mn.add(mi);

	}

	private void poll()
	{
		{
			if (timer == null)
				polltimer = new PollTimer(xw);
			polltimer.start();
		}
	}

	protected void addCommand(MenuItem mi, String command, final IJRequirement... requirements)
	{
		mi.setActionCommand(command);
		for (IJRequirement requirement : requirements)
		{
			if (requirement == IJRequirement.STACK && !xw.isStack())
				mi.setEnabled(false);
			if (requirement == IJRequirement.VOLUME && !xw.isVolume())
				mi.setEnabled(false);
		}
		mi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{

				try
				{
					MenuItem mi = (MenuItem) e.getSource();
					
					String command = mi.getActionCommand();
					runCommand(command, requirements);
				}
				catch (Exception ex)
				{
					XmippDialog.showInfo(null, ex.getMessage());
				}
			}
		});
	}//function addCommand
	
	/** Run ImageJ command */
	public void runCommand(String command, IJRequirement[] requirements){
		if (requirements != null)
			for (IJRequirement requirement : requirements)
				switch (requirement)
				{
				case IMAGEJ:
					XmippUtil.showImageJ(Tool.VIEWER);
					break;
				case BINARY:
					IJ.run("Make Binary");
					XmippDialog.showInfo(null, "make binary applied");
					break;
				case EIGHTBIT:
					IJ.run("8-bit");
					XmippDialog.showInfo(null, "8-bit applied");
					break;
				case THIRTYTWOBIT:
					IJ.run("32-bit");
					XmippDialog.showInfo(null, "32-bit applied");
					break;
				case RGB:
					IJ.run("RGB Color");
					XmippDialog.showInfo(null, "RGB color applied");
					break;
				
				}
		IJ.run(xw.getImagePlusLoader().getImagePlus(), command, "");
	}//function runCommand
	
	public Map<String, String> getKeyAssist()
	{
		Map<String, String> map = Collections.synchronizedMap(new LinkedHashMap<String, String>());
		map.put("Shift + Scroll Up", "Zoom in");
		map.put("Shift + Scroll Down", "Zoom out");
		map.put("Right click + Mouse move", "Moves image previously expanded");
		return map;
	}

}
