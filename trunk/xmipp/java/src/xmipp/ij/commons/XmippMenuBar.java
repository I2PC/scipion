/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package xmipp.ij.commons;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;

import java.awt.CheckboxMenuItem;
import java.awt.Menu;
import java.awt.MenuBar;
import java.awt.MenuItem;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;

import xmipp.jni.Filename;
import xmipp.utils.XmippFileChooser;
import javax.swing.JOptionPane;

import xmipp.viewer.windows.ImagesWindowFactory;

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
	private MenuItem openwithvv3ds;
	private MenuItem openwithvolumej;
	private Menu profilemn;
	private CheckboxMenuItem pollmi;
	protected Object timer;
	private final XmippIJWindow xw;
	private PollTimer polltimer;
	private MenuItem refreshmi;

	enum IJRequirement
	{
		BINARY, EIGHTBIT, IMAGEJ, STACK, THIRTYTWOBIT, RGB

	};

	public XmippMenuBar(XmippIJWindow xw)
	{
		this.xw = xw;
		// menubar menus
		filemn = new Menu("File");
		imagemn = new Menu("Image");
		advancedmn = new Menu("Advanced");

		add(filemn);
		add(imagemn);
		add(advancedmn);

		// menubar file menu
		savemi = new MenuItem("Save");
		savemi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				try
				{
					XmippMenuBar.this.xw.saveData();
				}
				catch (Exception ex)
				{
					JOptionPane.showMessageDialog(null, ex.getMessage());
				}

			}
		});
		saveasmi = new MenuItem("Save As...");
		saveasmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
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
					JOptionPane.showMessageDialog(null, ex.getMessage());
				}

			}
		});

		openwith3dmi = new MenuItem("Open with 3D Viewer");
		
//		openwith3dmi.setEnabled(Filename.isVolume(xw.getImageFilePath()));
		openwith3dmi.addActionListener(new ActionListener()
		{
			
			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				ImagePlus imp = XmippMenuBar.this.xw.getImagePlus();
				if (imp.getImageStackSize() == 1)
					JOptionPane.showMessageDialog(null, "Only for Stack");
				else
					ImagesWindowFactory.openImagePlusAs3D(imp);
				
			}
		});
		filemn.add(openwith3dmi);
		addIJMenuItem(filemn, "Open with Volume Viewer/3D Slicer", "Volume Viewer");
		addIJMenuItem(filemn, "Open with VolumeJ", "VolumeJ ");
		refreshmi = new MenuItem("Refresh");
		refreshmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				XmippMenuBar.this.xw.loadData();

			}
		});

		pollmi = new CheckboxMenuItem("Poll");
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

		filemn.add(savemi);
		filemn.add(saveasmi);
		addIJMenuItem(filemn, "Duplicate", "Duplicate...", IJRequirement.IMAGEJ);
		filemn.add(refreshmi);
		filemn.add(pollmi);

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

		addIJMenuItem(adjustmn, "Reslice", "Reslice [/]...");

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
				XmippIJUtil.showImageJ(Tool.VIEWER);
			}
		});

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
		mi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{

				try
				{
					String command = ((MenuItem) e.getSource()).getActionCommand();
					if (requirements != null)
						for (IJRequirement requirement : requirements)
							switch (requirement)
							{
							case IMAGEJ:
								XmippIJUtil.showImageJ(Tool.VIEWER);
								break;
							case BINARY:
								IJ.run("Make Binary");
								JOptionPane.showMessageDialog(null, "make binary applied");
								break;
							case EIGHTBIT:
								IJ.run("8-bit");
								JOptionPane.showMessageDialog(null, "8-bit applied");
								break;
							case THIRTYTWOBIT:
								IJ.run("32-bit");
								JOptionPane.showMessageDialog(null, "32-bit applied");
								break;
							case RGB:
								IJ.run("RGB Color");
								JOptionPane.showMessageDialog(null, "RGB color applied");
								break;
							case STACK:
								if (WindowManager.getCurrentImage().getImageStackSize() == 1)
									JOptionPane.showMessageDialog(null, "Only for Stack");
								return;
							}
					IJ.run(xw.getImagePlus(), command, "");
				}
				catch (Exception ex)
				{
					JOptionPane.showMessageDialog(null, ex.getMessage());
				}
			}
		});
	}

}
