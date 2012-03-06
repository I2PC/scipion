/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package xmipp.ij;

import ij.IJ;

import java.awt.CheckboxMenuItem;
import java.awt.Menu;
import java.awt.MenuBar;
import java.awt.MenuItem;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Timer;
import java.util.TimerTask;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import xmipp.particlepicker.ParticlePickerJFrame;

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
	private MenuItem openwithvv3ds;
	private MenuItem openwithvolumej;
	private MenuItem maskmi;
	private MenuItem thresholdmi;
	private MenuItem otsuthresholdmi;
	private MenuItem multiotsuthresholdmi;
	private MenuItem maxentropythresholdmi;
	private MenuItem mixturemodthresholdmi;
	private MenuItem voxelcountermi;
	private MenuItem erodemi;
	private MenuItem dilatemi;
	private MenuItem openmi;
	private MenuItem closemi;
	private MenuItem floatmorphomi;
	private MenuItem outlinemi;
	private MenuItem skeletonizemi;
	private MenuItem fillholesmi;
	private MenuItem distancemapmi;
	private MenuItem ultimatepointsmi;
	private MenuItem watershedmi;
	private MenuItem voronoimi;
	private Menu profilemn;
	private MenuItem lineanalyzermi;
	private MenuItem ovalpplotmi;
	private MenuItem radialpplotanglemi;
	private MenuItem radialpplotheightmi;
	private MenuItem contourplottermi;
	private MenuItem brightcontrastmi;
	private MenuItem enhancecontrastmi;
	private MenuItem substractbgmi;
	private MenuItem gaussianblurmi;
	private MenuItem convolvemi;
	private MenuItem medianmi;
	private MenuItem fftmi;
	private CheckboxMenuItem pollmi;
	protected Object timer;
	private final XmippIJWindow xw;
	private PollTimer polltimer;
	private MenuItem refreshmi;

	public XmippMenuBar(XmippIJWindow xw)
	{
		this.xw = xw;
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
				JFileChooser fc = new JFileChooser();
				int returnVal = fc.showOpenDialog(null);

				try
				{
					if (returnVal == JFileChooser.APPROVE_OPTION)
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
		openwithvv3ds = new MenuItem("Open with Volume Viewer/3D Slicer");
		openwithvolumej = new MenuItem("Open with VolumeJ");
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
		filemn.add(openwith3dmi);
		filemn.add(openwithvv3ds);
		filemn.add(openwithvolumej);
		filemn.add(refreshmi);
		filemn.add(pollmi);
		

		// menubar image menu
		infomn = new Menu("Info");
		adjustmn = new Menu("Adjust");
		transformmn = new Menu("Transform");
		filtersmn = new Menu("Filters");
		maskmi = new MenuItem("Mask Tool Bar");
		addCommand(maskmi, "Mask Tool Bar");

		imagemn.add(infomn);
		imagemn.add(adjustmn);
		imagemn.add(transformmn);
		imagemn.add(filtersmn);
		imagemn.add(maskmi);

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
		addCommand(bandpassmi, "Bandpass Filter...");
		admi = new MenuItem("Anisotropic Diffusion");
		addCommand(admi, "Anisotropic Diffusion...");
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
		profilemn = new Menu("Profile");

		advancedmn.add(imagejmi);
		advancedmn.add(thresholdingmn);
		advancedmn.add(binarymn);
		advancedmn.add(processmn);
		advancedmn.add(drawmn);
		advancedmn.add(profilemn);

		// advanced threshold menu
		thresholdmi = new MenuItem("Threshold");
		otsuthresholdmi = new MenuItem("Otsu Threshold");
		multiotsuthresholdmi = new MenuItem("Multi Otsu Threshold");
		maxentropythresholdmi = new MenuItem("Maximum Entropy Threshold");
		mixturemodthresholdmi = new MenuItem("Mixture Modeling Threshold");

		thresholdingmn.add(thresholdmi);
		thresholdingmn.add(otsuthresholdmi);
		thresholdingmn.add(multiotsuthresholdmi);
		thresholdingmn.add(maxentropythresholdmi);
		thresholdingmn.add(mixturemodthresholdmi);

		// advanced binary menu
		voxelcountermi = new MenuItem("Voxel Counter");
		addCommand(voxelcountermi, "Voxel Counter");
		erodemi = new MenuItem("Erode");
		addCommand(erodemi, "Erode");
		dilatemi = new MenuItem("Dilate");
		openmi = new MenuItem("Open");
		closemi = new MenuItem("Close");
		floatmorphomi = new MenuItem("Float Morphology");
		outlinemi = new MenuItem("Outline");
		fillholesmi = new MenuItem("Fill Holes");
		skeletonizemi = new MenuItem("Skeletonize");
		distancemapmi = new MenuItem("Distance Map");
		ultimatepointsmi = new MenuItem("Ultimate Points");
		watershedmi = new MenuItem("Water Shed");
		voronoimi = new MenuItem("Voronoi");

		binarymn.add(voxelcountermi);
		binarymn.add(erodemi);
		binarymn.add(dilatemi);
		binarymn.add(openmi);
		binarymn.add(closemi);
		binarymn.add(floatmorphomi);
		binarymn.add(outlinemi);
		binarymn.add(fillholesmi);
		binarymn.add(skeletonizemi);
		binarymn.add(watershedmi);
		binarymn.add(voronoimi);

		// advanced process menu
		brightcontrastmi = new MenuItem("Brightness/Contrast");
		enhancecontrastmi = new MenuItem("Enhance Contrast");
		substractbgmi = new MenuItem("Substract Background");
		gaussianblurmi = new MenuItem("Gaussian Blur");
		convolvemi = new MenuItem("Convolve");
		medianmi = new MenuItem("Median");
		fftmi = new MenuItem("FFT");

		processmn.add(brightcontrastmi);
		processmn.add(enhancecontrastmi);
		processmn.add(substractbgmi);
		processmn.add(gaussianblurmi);
		processmn.add(convolvemi);
		processmn.add(medianmi);
		processmn.add(fftmi);

		// advanced drawn menu

		// advanced profile menu
		lineanalyzermi = new MenuItem("Line Analyzer");
		ovalpplotmi = new MenuItem("Oval Profile Plot");
		radialpplotanglemi = new MenuItem("Radial Profile Plot Angle");
		radialpplotheightmi = new MenuItem("Radial Profile Plot Height");
		contourplottermi = new MenuItem("Contour Plotter");

		profilemn.add(lineanalyzermi);
		profilemn.add(ovalpplotmi);
		profilemn.add(radialpplotanglemi);
		profilemn.add(radialpplotheightmi);
		profilemn.add(contourplottermi);

		imagejmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				XmippImageWindow.openImageJ(Tool.VIEWER);
			}
		});

	}
	


	private void poll()
	{
		{
			if(timer == null)
				polltimer = new PollTimer(xw);
			polltimer.start();
		}
	}

	protected void addCommand(MenuItem mi, String command)
	{
		addCommand(mi, command, false);
	}

	protected void addCommand(MenuItem mi, String command, boolean isrequiredij)
	{
		if (isrequiredij)
			requireij.add(command);
		mi.setActionCommand(command);
		mi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{

				String command = ((MenuItem) e.getSource()).getActionCommand();
				if (requireij.contains(command))
					XmippImageWindow.openImageJ(Tool.VIEWER);
				if (command.equals("Anisotropic Diffusion..."))
					IJ.run("8-bit");
				IJ.run(command);
			}
		});
	}

}
