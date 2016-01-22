/*
cd * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package xmipp.ij.commons;

import ij.CommandListener;
import ij.Executer;
import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.plugin.frame.Recorder;
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
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Color3f;

import xmipp.utils.QuickHelpJDialog;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippWindowUtil;

/**
 * 
 * @author Juanjo Vega
 */
public class XmippMenuBar extends MenuBar
{

	protected final XmippIJWindow xw;
	protected Menu filemn;
	protected Menu imagemn;
	protected Menu advancedmn;
	protected MenuItem savemi;
	protected MenuItem saveasmi;
	protected MenuItem openwith3dmi;
	protected Menu infomn;
//	protected Menu adjustmn;
	protected Menu transformmn;
//	protected Menu filtersmn;
	protected Menu thresholdingmn;
	protected Menu binarymn;
	protected Menu processmn;
	protected Menu drawmn;
	protected MenuItem imagejmi;
	protected Menu profilemn;
	protected CheckboxMenuItem pollmi;
	protected Object timer;
	protected PollTimer polltimer;
	protected MenuItem refreshmi;
	protected CheckboxMenuItem wrapmi;
	protected CheckboxMenuItem ugmi;
	protected MenuItem exitmi;
	protected Menu helpmn;
	protected MenuItem onlinemi;
	protected MenuItem keyassistmi;
	protected QuickHelpJDialog keyassistdlg;
	protected boolean invertx;
	protected boolean inverty;
	protected CheckboxMenuItem activemi;
	protected String command;
	protected List<IJCommand> filters;
	protected ImagePlusLoader iploader;

	enum IJRequirement
	{
		BINARY, EIGHTBIT, IMAGEJ, STACK, THIRTYTWOBIT, RGB, VOLUME, INVERTX, INVERTY

	};

	public XmippMenuBar(XmippIJWindow xw)
	{
		filters = new ArrayList<IJCommand>();
		this.xw = xw;
		this.iploader = xw.getImagePlusLoader();
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
//		savemi = new MenuItem("Save");
//		savemi.setShortcut(new MenuShortcut(KeyEvent.VK_S));
//		savemi.addActionListener(new ActionListener()
//		{
//
//			@Override
//			public void actionPerformed(ActionEvent e)
//			{
//				try
//				{
//					if(XmippMenuBar.this.iploader.existsFile())
//						XmippMenuBar.this.xw.saveData();
//					else
//						XmippMenuBar.this.saveAs();
//				}
//				catch (Exception ex)
//				{
//					ex.printStackTrace();
//					XmippDialog.showInfo(null, ex.getMessage());
//				}
//
//			}
//		});
		saveasmi = new MenuItem("Save As...");
		saveasmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				saveAs();
			}
		});

//		openwith3dmi = new MenuItem("Open with 3D Viewer");
//		openwith3dmi.setEnabled(xw.isVolume());
//		openwith3dmi.setEnabled(Filename.isVolume(xw.getImageFilePath()));
//		openwith3dmi.addActionListener(new ActionListener()
//		{
//			
//			@Override
//			public void actionPerformed(ActionEvent arg0)
//			{
//				ImagePlus imp = XmippMenuBar.this.iploader.getImagePlus();
//				if (imp.getImageStackSize() == 1)
//					XmippDialog.showInfo(null, "Only for Stack");
//				else
//					openImagePlusAs3D(imp);
//				
//			}
//		});
//		filemn.add(openwith3dmi);
		MenuItem volviewermi = new MenuItem("Open with Volume Viewer/3D Slicer");
		volviewermi.setActionCommand("Volume Viewer");
		volviewermi.setEnabled(xw.isVolume());
		volviewermi.addActionListener(new ActionListener()
		{
			
			@Override
			public void actionPerformed(ActionEvent e)
			{
				MenuItem mi = (MenuItem) e.getSource();
				runCommand(mi.getActionCommand());
				
			}
		});
		filemn.add(volviewermi);
//		addIJMenuItem(filemn, "Open with Volume Viewer/3D Slicer", "Volume Viewer", IJRequirement.VOLUME);
//		addIJMenuItem(filemn, "Open with VolumeJ", "VolumeJ ", IJRequirement.VOLUME);
		refreshmi = new MenuItem("Refresh");
		refreshmi.setEnabled(iploader.allowsPoll());
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
		pollmi.setEnabled(iploader.allowsPoll());
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
		ugmi.setEnabled(iploader.hasGeometry());
		ugmi.setState(iploader.getUseGeometry());
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
		wrapmi.setEnabled(iploader.hasGeometry());
		wrapmi.setState(iploader.isWrap());
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



//		filemn.add(savemi);
		filemn.add(saveasmi);
//		addIJMenuItem(filemn, "Duplicate", "Duplicate...", IJRequirement.IMAGEJ);
		filemn.add(refreshmi);
		filemn.add(pollmi);
		filemn.add(ugmi);
		filemn.add(wrapmi);
		filemn.add(exitmi);

		// menubar image menu
//		infomn = new Menu("Info");
//		adjustmn = new Menu("Adjust");
//		filtersmn = new Menu("Filters");

//		imagemn.add(infomn);
//		imagemn.add(adjustmn);
//		imagemn.add(transformmn);
//		imagemn.add(filtersmn);
		

		// image info menu
//		addIJMenuItem(infomn, "Show Info", "Show Info...");
//		addIJMenuItem(infomn, "Properties", "Properties...");
//		addIJMenuItem(infomn, "Histogram", "Histogram", IJRequirement.IMAGEJ);
//		addIJMenuItem(infomn, "Plot Profile", "Plot Profile", IJRequirement.IMAGEJ);

		// image adjust menu

//		addIJMenuItem(adjustmn, "Crop", "Crop", IJRequirement.IMAGEJ);
//		addIJMenuItem(adjustmn, "Scale", "Scale...");
//		addIJMenuItem(adjustmn, "Untilt Stack", "Untilt Stack", IJRequirement.STACK);

//		addIJMenuItem(adjustmn, "Reslice", "Reslice [/]...", IJRequirement.VOLUME);

		// image transform menu
		
//		addIJMenuItem(transformmn, "Rotate 90 Degrees Left", "Rotate 90 Degrees Left");
//		addIJMenuItem(transformmn, "Rotate 90 Degrees Right", "Rotate 90 Degrees Right");

		// image filters menu
		CheckboxMenuItem gbmi = addIJMenuItem(imagemn, XmippImageJ.gaussianBlurFilter);
		addIJMenuItem(imagemn, XmippImageJ.bandPassFilter);
//		addIJMenuItem(imagemn, XmippImageJ.anisotropicDiffFilter, IJRequirement.EIGHTBIT);
		CheckboxMenuItem ecmi = addIJMenuItem(imagemn, XmippImageJ.enhanceContrastFilter);
		addIJMenuItem(imagemn, XmippImageJ.brightnessContrastFilter);
		addIJMenuItem(imagemn, XmippImageJ.invertLUTFilter);
		addIJMenuItem(imagemn, XmippImageJ.substractBackgroundFilter);
		// menubar advanced menu
		imagejmi = new MenuItem("ImageJ");
		imagejmi.setShortcut(new MenuShortcut(KeyEvent.VK_I));
		transformmn = new Menu("Transform");
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
//		addIJMenuItem(advancedmn, "Masks Tool Bar", "Masks Tool Bar", IJRequirement.IMAGEJ);// missing
		// plugin
		onlinemi = new MenuItem("Online Help");
		onlinemi.addActionListener(new ActionListener()
		{
			
			

			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				XmippWindowUtil.openURI("https://github.com/biocompwebs/scipion/wiki/ShowJ");
			}
		});
		helpmn.add(onlinemi);
		keyassistmi = new MenuItem("Tips");
		keyassistmi.addActionListener(new ActionListener()
		{
			
			

			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				try
				{
				if(keyassistdlg == null)
					keyassistdlg = new QuickHelpJDialog(null, false, "Tips...", getKeyAssist());
				keyassistdlg.setVisible(true);
				}
				catch(Exception e)
				{
					XmippDialog.showInfo(null, e.getMessage());
				}
			}
		});
		helpmn.add(keyassistmi);
		
		addIJMenuItem(transformmn, "Flip Horizontally", "Flip Horizontally", IJRequirement.INVERTX);
		addIJMenuItem(transformmn, "Flip Vertically", "Flip Vertically", IJRequirement.INVERTY);

		// advanced threshold menu
		addIJMenuItem(thresholdingmn, "Threshold", "Threshold...");
		addIJMenuItem(thresholdingmn, "Otsu Threshold", "Otsu Thresholding", IJRequirement.EIGHTBIT);
		addIJMenuItem(thresholdingmn, "Multi Otsu Threshold", "Multi OtsuThreshold");
		addIJMenuItem(thresholdingmn, "Maximum Entropy Threshold", "Entropy Threshold", IJRequirement.EIGHTBIT);
		addIJMenuItem(thresholdingmn, "Mixture Modeling Threshold", "Mixture Modeling", IJRequirement.EIGHTBIT);
		addIJMenuItem(thresholdingmn, "Robust Automatic Threshold Selection", "RATS ");
		//addIJMenuItem(thresholdingmn, "Simple Iterative Object Extraction", "SIOX Segmentation", IJRequirement.RGB);

		// advanced binary menu
		addIJMenuItem(binarymn, "Voxel Counter", "Voxel Counter", IJRequirement.STACK);// stack required
		addIJMenuItem(binarymn, "Erode", "Erode", IJRequirement.BINARY);
		addIJMenuItem(binarymn, "Dilate", "Dilate", IJRequirement.BINARY);
		addIJMenuItem(binarymn, "Open", "Open", IJRequirement.BINARY);
		addIJMenuItem(binarymn, "Close", "Close", IJRequirement.BINARY);
		addIJMenuItem(binarymn, "Float Morphology", "Float Morphology", IJRequirement.THIRTYTWOBIT);// missing plugin
		addIJMenuItem(binarymn, "Outline", "Outline", IJRequirement.BINARY, IJRequirement.EIGHTBIT);
		addIJMenuItem(binarymn, "Fill Holes", "Fill Holes", IJRequirement.BINARY, IJRequirement.EIGHTBIT);
		addIJMenuItem(binarymn, "Skeletonize", "Skeletonize", IJRequirement.BINARY, IJRequirement.EIGHTBIT);
		addIJMenuItem(binarymn, "Distance Map", "Distance Map", IJRequirement.BINARY, IJRequirement.EIGHTBIT);
		addIJMenuItem(binarymn, "Ultimate Points", "Ultimate Points", IJRequirement.BINARY, IJRequirement.EIGHTBIT);
		addIJMenuItem(binarymn, "Watershed", "Watershed");
		addIJMenuItem(binarymn, "Voronoi", "Voronoi", IJRequirement.BINARY, IJRequirement.EIGHTBIT);

		// advanced process menu

		
		addIJMenuItem(processmn, "Convolve", "Convolve...");
		addIJMenuItem(processmn, "Median", "Median...");
		addIJMenuItem(processmn, "FFT", "FFT", IJRequirement.IMAGEJ);// memory error
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
		addFilterAppliedListener();
		ImagePlus imp = iploader.getImagePlus();
		if(xw.getParams() != null && imp.getWidth() > 1024 && imp.getStackSize() == 1)
		{
            String value = System.getenv("SCIPION_MIC_NOGAUSS");
                        
            if (value == null || !value.equals("1"))
            {
                gbmi.setState(true);
                IJCommand ijcmd = new IJCommand(XmippImageJ.gaussianBlurFilter, "sigma =2");
                IJ.run(imp, ijcmd.getCommand(), ijcmd.getOptions());
                filters.add(ijcmd);
            }
            
            value = System.getenv("SCIPION_MIC_NOENHANCE");            
            if (value == null || !value.equals("1"))
            {
                ecmi.setState(true);
                IJCommand ijcmd = new IJCommand(XmippImageJ.enhanceContrastFilter, "saturated=0.4");
                IJ.run(imp, ijcmd.getCommand(), ijcmd.getOptions());
                filters.add(ijcmd);
            }
		}
			
	}
	
	protected void addFilterAppliedListener() {

        Recorder.record = true;
     // detecting if a command is thrown by ImageJ
        Executer.addCommandListener(new CommandListener() {
            public String commandExecuting(String command) {


                XmippMenuBar.this.command = command;
                return command;

            }
        });

        // detecting if a command is thrown by ImageJ
        
        ImagePlus.addImageListener(new ImageListener() {

            @Override
            public void imageUpdated(ImagePlus imp) {
            	//activemi will contain last filter selected, if its the command runned we confirm selection
            	if(activemi != null && activemi.getActionCommand().equals(command))
            	{
                   activemi.setState(true);
                   String options = Recorder.getCommandOptions();
                   filters.add(new IJCommand(command, options));
                   activemi = null;
            	}
                
            }

            @Override
            public void imageOpened(ImagePlus arg0) {

            }

            @Override
            public void imageClosed(ImagePlus arg0) {
                // TODO Auto-generated method stub

            }
        });
	}
        

    protected void close()
    {
        Frame w = (Frame)xw;
        w.setVisible(false);
        w.dispose();
        XmippApplication.removeInstance(false);
    }

	
	protected void useGeometry(boolean ug)
	{
		xw.getImagePlusLoader().setUseGeometry(ug);
		xw.loadData();
	}

	protected void saveAs()
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
	
	protected CheckboxMenuItem addIJMenuItem(Menu mn, String command, IJRequirement... requirements)
	{
		return addIJMenuItem(mn, command, command, requirements);

	}

	protected CheckboxMenuItem addIJMenuItem(Menu mn, String name, String command, IJRequirement... requirements)
	{
		CheckboxMenuItem mi = new CheckboxMenuItem(name);
		addCommand(mi, command, requirements);
		mn.add(mi);
		return mi;
	}

	protected void poll()
	{
			if (timer == null)
				polltimer = new PollTimer(xw);
			polltimer.start();
	}

	protected void addCommand(CheckboxMenuItem mi, String command, final IJRequirement... requirements)
	{
		mi.setActionCommand(command);
		for (IJRequirement requirement : requirements)
		{
			if (requirement == IJRequirement.STACK && !xw.isStack())
				mi.setEnabled(false);
			if (requirement == IJRequirement.VOLUME && !xw.isVolume())
				mi.setEnabled(false);
		}
		mi.addItemListener(new ItemListener()
		{

			@Override
			public void itemStateChanged(ItemEvent e)
			{
				try
				{
					CheckboxMenuItem mi = (CheckboxMenuItem) e.getSource();
					String command = mi.getActionCommand();
					if(mi.getState())//checked
					{
						activemi = mi;
						mi.setState(false);//confirm selected only after filter applied
						runCommand(command, requirements);
						
					}
					else
					{
						filters.remove(getFilter(command));
						xw.loadData();
						applyFilters();
					}
					
				}
				catch (Exception ex)
				{
					ex.printStackTrace();
				}
			}
		});
	}//function addCommand
	
	public void applyFilters()
	{
		ImagePlus imp = xw.getImagePlusLoader().getImagePlus();
		for(IJCommand filter: filters)
			IJ.run(imp, filter.getCommand(), filter.getOptions());
	}
	
	
	public IJCommand getFilter(String command)
	{
		for(IJCommand filter: filters)
		{
			if(filter.getCommand().equals(command))
				return filter;
		}
		return null;
	}
	/** Run ImageJ command */
	public void runCommand(String command, IJRequirement... requirements){
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
				case INVERTX:
					invertx = !invertx;
					xw.getCanvas().setInvertX(invertx);
					break;
				case INVERTY:
					inverty = !inverty;
					xw.getCanvas().setInvertY(inverty);
					break;
				
				}
		ImagePlus imp = xw.getImagePlusLoader().getImagePlus();
		IJ.run(imp, command, "");
	}//function runCommand
	
	public Map<Object, Object> getKeyAssist()
	{
		Map<Object, Object> map = Collections.synchronizedMap(new LinkedHashMap<Object, Object>());
		map.put("Shift + scroll up", "Zoom in");
		map.put("Shift + scroll down", "Zoom out");
		map.put("Right click + mouse move", "Moves image previously expanded");
		return map;
	}
        
        
}
