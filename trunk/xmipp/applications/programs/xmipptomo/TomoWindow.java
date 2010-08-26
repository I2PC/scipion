/***************************************************************************
 *
 * @author: Jesus Cuenca (jcuenca@cnb.csic.es)
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

package xmipptomo;

import ij.IJ;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.ImageLayout;
import ij.io.OpenDialog;
import ij.io.SaveDialog;

import java.awt.*;
import java.awt.event.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.FileNotFoundException;
import java.io.IOException;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;


/**
 * @author jcuenca
 * Custom visualization window, using Swing and MVC paradigm
 * 
 * implements: WindowListener, AdjustmentListener,MouseMotionListener,ActionListener, PropertyChangeListener to handle GUI events
 * 			DialogListener to capture GenericDialogs
 */
/** Protocol for adding new buttons to the workflow/UI:
 *  - set the label in Buttons enum as static String, and its corresponding ImageJ command (if any)
 *  - run addButton 
 *  - update actionPerformed method
 */
public class TomoWindow extends JFrame implements WindowListener, AdjustmentListener,MouseMotionListener,ActionListener, PropertyChangeListener,DialogListener{
	// for serialization only - just in case
	private static final long serialVersionUID = -4063711975454855701L;
	
	// total number of  digits to display as a String
	private static int MAX_DIGITS = 9;
	// Window title (without extra data; see getTitle() )
	private final static String TITLE="XmippTomo";
	// minimum dimension of tabbed menu panel (window top)
	private static int MENUPANEL_MINWIDTH=400, MENUPANEL_MINHEIGHT=100;
	// miliseconds to wait for a Plugin dialog to display (so the dialog capture can start)
	private static int WAIT_FOR_DIALOG=800;
	
	// labels of buttons and their associated commands (if any)
	private static enum Buttons {
		CMD_LOAD("Load",""),CMD_SAVE("Save",""),
		CMD_GAUSSIAN("Gaussian",GaussianPlugin.COMMAND),CMD_MEDIAN("Median", MedianPlugin.COMMAND),CMD_SUB_BACKGROUND("Bandpass - substract background",LowPassPlugin.COMMAND),
		CMD_APPLY("Apply",""),
		CMD_PRINT_WORKFLOW("Print workflow","");

		private final String label;
		private final String imageJCmd;
		
		Buttons(String l,String cmd) { this.label=l;this.imageJCmd=cmd;}
		public String label(){return label;}
		public String imageJCmd(){return imageJCmd;}
	};
	
	// Menu labels
	private enum Menu {
		FILE("File"),PREPROC("Preprocess"),	ALIGN("Align"), DEBUG("Debug");

		private String label;

		private Menu(String l) { this.label=l;}
		public final String toString(){return label;}
	}
	
	private boolean closed=true;
	// true while reloading a file that was resized
	private boolean reloadingFile=false;
	// true if changes are saved (so you can close the window without pain)
	private boolean changeSaved=true;
	// window identifier - useful when you've got more than 1 window
	private int windowId=-1;
	
	/* Window GUI components */
	/* 4 main panels */
	private JTabbedPane menuPanel;
	private JPanel viewPanel,controlPanel,statusPanel;
	
	// Text fields and labels
	private JTextField tiltTextField;
	private JLabel statusLabel;
	
	// The canvas where current projection is displayed as an image
	private ImageCanvas ic;
	
	// Control scrollbars
	private LabelScrollbar projectionScrollbar;
	
	// the model stores the data that this window shows - @see TomoData
	private TomoData model=null;
	
	Point cursorLocation=new Point();
	
	// action history hints - so we can navigate the workflow
	private UserAction firstAction,lastAction;
	
	// current plugin - DialogCapture thread uses it (when dialotItemChanged is called)
	private Plugin plugin;
	
	// ImportDataThread: load & display projections in parallel
	private class ImportDataThread extends Thread{
		private TomoData dataModel;
		private boolean resize=false;
		
		ImportDataThread(TomoData model,boolean resize){
			dataModel=model;
			this.resize=resize;
		}
		
		public void run() {
			try{
				new TiltSeriesOpener(resize).read(dataModel);
			}catch (FileNotFoundException ex){
				Xmipp_Tomo.debug("ImportDataThread.run - " + ex.getMessage());
			}catch (IOException ex){
				Xmipp_Tomo.debug("ImportDataThread.run - Error opening file");
			}catch (InterruptedException ex){
				Xmipp_Tomo.debug("ImportDataThread.run - Interrupted exception");
			}catch (Exception ex){
				Xmipp_Tomo.debug("ImportDataThread.run - unexpected exception", ex);
			}
		}
	}
	
	// implement UI concurrency with private class that extends Thread
	// ActionThread: handle user actions in parallel (no GUI locking)
	private class ActionThread extends Thread{
		private TomoWindow tw;
		private String command;
		
		ActionThread(TomoWindow window,String command){
			tw=window;
			this.command=command;
		}
		
		public void run() {
			// select proper action method based on button's label
			if (command.equals(Buttons.CMD_LOAD.label())){
				tw.actionLoad();
			}else if (command.equals(Buttons.CMD_GAUSSIAN.label())){
				tw.actionRunIjCmd(Buttons.CMD_GAUSSIAN.label(), Buttons.CMD_GAUSSIAN.imageJCmd());
			}else if (command.equals(Buttons.CMD_MEDIAN.label())){
				tw.actionRunIjCmd(Buttons.CMD_MEDIAN.label(), Buttons.CMD_MEDIAN.imageJCmd());
			}else if (command.equals(Buttons.CMD_SUB_BACKGROUND.label())){
				tw.actionRunIjCmd(Buttons.CMD_SUB_BACKGROUND.label(), Buttons.CMD_SUB_BACKGROUND.imageJCmd());
			}else if (command.equals(Buttons.CMD_APPLY.label())){
				tw.actionApply();
			}else if (command.equals(Buttons.CMD_PRINT_WORKFLOW.label())){
				Xmipp_Tomo.printWorkflow();
			}
		}
	}
	
	// implement UI concurrency with private class that extends Thread
	// DialogCaptureThread: intercept generic dialogs
	private class DialogCaptureThread extends Thread{
		private TomoWindow window;
		
		DialogCaptureThread(TomoWindow w){
			window=w;
		}
		
		public void run() {
			try{
				sleep(WAIT_FOR_DIALOG);
				window.captureDialog();
			}catch (Exception ex){
				Xmipp_Tomo.debug("DialogCaptureThread.run - unexpected exception", ex);
			}
		}
	}
	
	/* METHODS ------------------------------------------------------------- */
	
	/* For the buggy GridBagLayout...
	// row: 1..N
	private void addToGrid(Container c,int row, int constraints_fill){
		GridBagConstraints constraints=new GridBagConstraints();
		constraints.gridx=0; constraints.gridy = row-1;
		constraints.anchor = GridBagConstraints.SOUTH;
		constraints.fill = constraints_fill;	
		//gbl.setConstraints(c, constraints);
		getContentPane().add(c,constraints); 
	}*/
	
	/**
	 * Create a window with its main panels (but no data/model available) and show it
	 */
	public TomoWindow(){
		// set this window as listener of general keyboard&mouse events
 		addWindowListener(this);
 		// EXIT_ON_CLOSE finishes ImageJ too...
 		// DO NOTHING allows for closing confirmation dialogs and the like
	    setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		WindowManager.addWindow(this);
		addMainPanels();	
	}
	
	public TomoWindow(int windowId){
		this();
		setWindowId(windowId);		
	}
	
	/** Add button to a panel
	 * @param label
	 */
	void addButton(String label,Container panel) {
		Button b = new Button(label);
		b.addActionListener(this);
		// b.addKeyListener(IJ.getInstance());
		panel.add(b); 
	}
	
	/** 
	 * Add the main panels with their components, and associate each text field with model's documents (if a model is already set)
	 */
	private void addMainPanels(){
		// Window title
 		setTitle(getTitle());
 		
 		setLayout(new BoxLayout(getContentPane(),BoxLayout.PAGE_AXIS));

		// BoxLayout component order is sequential, so all the panels
 		// MUST be added in order (even if they are empty)
 		
		// TABBED MENU PANEL
		menuPanel= new JTabbedPane();
		addMenuTabs();	
		menuPanel.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
		menuPanel.setPreferredSize(new Dimension(MENUPANEL_MINWIDTH,MENUPANEL_MINHEIGHT));
		add(menuPanel);
		
		// VIEW & CONTROLS PANELS
		viewPanel = new JPanel();
		add(viewPanel);
		
		controlPanel = new JPanel();
		add(controlPanel);
		
		if(getModel() != null){
			addView();
			addControls();
		}

		// STATUS PANEL
		statusPanel = new JPanel();
		FlowLayout fl=new FlowLayout();
		fl.setAlignment(FlowLayout.LEFT);
		statusPanel.setLayout(fl);
		statusPanel.setMaximumSize(new Dimension(Short.MAX_VALUE,50));

		statusLabel=new JLabel("Ready");
		statusPanel.add(statusLabel);
		updateStatusText();
		add(statusPanel);
		
		closed=false;
		pack(); // adjust GUI size
	}
	
	/**
	 * It requires that the model has been set (with setModel() ), and that at least 1 image has been loaded
	 * @param model
	 */
	private void addView(){
		if(getModel() == null)
			return;
		
		//remove previous canvas if any
		for(int i=0; i < viewPanel.getComponentCount(); i++)
			viewPanel.remove(i);
		
		ic=new ImageCanvas(getModel().getImage());
		viewPanel.setLayout(new ImageLayout(ic));
		viewPanel.add(ic);
		ic.addMouseMotionListener(this);
	
		pack(); // adjust window size to host the new view panel
	}
	
	/**
	 * It also requires that the model has been set (with setModel() ), and that basic info (number of projections,...) has been loaded
	 * @param model
	 */
	private void addControls(){
		if(getModel() == null)
			return;
			
		//remove previous controls if any
		for(int i=0; i < controlPanel.getComponentCount(); i++)
			controlPanel.remove(i);
		
		FlowLayout fl2=new FlowLayout();
		fl2.setAlignment(FlowLayout.CENTER);
		controlPanel.setLayout(fl2);

		
		JLabel tiltTextLabel=new JLabel("Tilt");
		tiltTextField=new JTextField(5);
		// the user won't change angles one by one - he will change all of them at once with the tilt change dialog
		tiltTextField.setEditable(false);
		controlPanel.add(tiltTextLabel);		
		controlPanel.add(tiltTextField);
		getModel().setTiltModel(tiltTextField.getDocument());
		
		projectionScrollbar= new LabelScrollbar(1,getModel().getNumberOfProjections());
		projectionScrollbar.setText("Projection #");
		projectionScrollbar.addAdjustmentListener(this);
		controlPanel.add(projectionScrollbar);
		
		pack(); // adjust window size to host the new controls panel
				
	}
	
	private void addMenuTabs(){
		// File
		// ButtonTabComponent file_tab=new ButtonTabComponent(menuPanel), preproc_tab=new ButtonTabComponent(menuPanel);
		JPanel fileTabPanel=new JPanel(false); // false = no double buffer (flicker)
		//menuPanel.setMnemonicAt(0, KeyEvent.VK_1);
		addButton(Buttons.CMD_LOAD.label(),fileTabPanel);
		addButton(Buttons.CMD_SAVE.label(),fileTabPanel);
		menuPanel.addTab(Menu.FILE.toString(),fileTabPanel); // try to add the tab when it's ready (all controls inside)

		// Preprocessing
		JPanel preprocTabPanel=new JPanel();
		addButton(Buttons.CMD_GAUSSIAN.label(),preprocTabPanel);
		addButton(Buttons.CMD_MEDIAN.label(),preprocTabPanel);
		addButton(Buttons.CMD_SUB_BACKGROUND.label(),preprocTabPanel);
		addButton(Buttons.CMD_APPLY.label(),preprocTabPanel);
		menuPanel.addTab(Menu.PREPROC.toString(),preprocTabPanel);
		
		// Debugging
		if(Xmipp_Tomo.TESTING == 1){
			JPanel debugTabPanel=new JPanel();
			addButton(Buttons.CMD_PRINT_WORKFLOW.label(),debugTabPanel);
			menuPanel.addTab(Menu.DEBUG.toString(),debugTabPanel);
		}
	}
	
	
	
	/**
	 * warning: naming this method "show" leads to infinite recursion (since setVisible calls show in turn)
	 */
	public void display(){
		 pack();
	     setVisible(true);
	}
	
	
	private void setStatus(String text){
		getStatusLabel().setText(text);
	}
	
	private void updateStatusText(){
		// current/total projections are shown in the projection scrollbar itself
		if(getModel() != null)
			setStatus("x = " + getCursorX() + ", y = " + getCursorY() + ", value = " + getCursorValueAsString());
	}
	
	int getCursorDistance(int x, int y){
		return (int) cursorLocation.distance(x,y);
	}
	
	private void refreshImageCanvas(){
		if(ic != null){
			ic.setImageUpdated();
			ic.repaint();
		}
	}
	
	/**
	 * Scrollbar events
	 * @param e
	 */
	public synchronized void adjustmentValueChanged(AdjustmentEvent e){
		// Xmipp_Tomo.debug("TomoWindow.adjustmentValueChanged"+e.getSource().toString());
		if(e.getSource()==projectionScrollbar){
			// scrollbar range is 1..N, like projections array
			getModel().setCurrentProjection(projectionScrollbar.getValue());
			refreshImageCanvas();
			updateStatusText(); // projection number
		}
		notify();
	}
	
	/* button/keyboard events - run in a different thread to not to block the GUI
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	public void actionPerformed(ActionEvent e) {
	
		String label = e.getActionCommand();
		if (label==null)
			return;
		
		new ActionThread(this, label).start();
		
	} // actionPerformed end
	
	/* -------------------------- Action management --------------------*/
	
	private void addUserAction(UserAction a){
		Xmipp_Tomo.addUserAction(getLastAction(), a);
		setLastAction(a);
	}
	
	private void actionLoad(){
		String path = dialogOpen();
		if((path == null) || ("".equals(path)))
			return;
		
		setModel(new TomoData(path));
	
		
		try{
			// import data in one thread and hold this thread until the first projection is loaded
			setReloadingFile(false);
			(new Thread(new ImportDataThread(getModel(),true))).start();
			getModel().waitForFirstImage();
		}catch (InterruptedException ex){
			Xmipp_Tomo.debug("Xmipp_Tomo - Interrupted exception");
		}

		setTitle(getTitle());
		addView();
		addControls();
		addUserAction(new UserAction(getWindowId(),"Load",getModel().getFileName()));
	}
	
	private void actionRunIjCmd(String label,String cmd){
		if(cmd == null)
			return;
		
		if(label != null)
			setStatus(label + " - started");
		
		setChangeSaved(false);
		
		if(cmd.equals(Buttons.CMD_GAUSSIAN.imageJCmd()))
			setPlugin(new GaussianPlugin());
		else if(cmd.equals(Buttons.CMD_MEDIAN.imageJCmd()))
			setPlugin(new MedianPlugin());
		else if(cmd.equals(Buttons.CMD_SUB_BACKGROUND.imageJCmd()))
			setPlugin(new LowPassPlugin());
		
		(new Thread(new DialogCaptureThread(this))).start();
		
		WindowManager.setTempCurrentImage(getModel().getImage());
		IJ.run(cmd);

		refreshImageCanvas();
		
		setStatus("Done");
		
		addUserAction(new UserAction(getWindowId(),cmd,getPlugin()));
		
		setPlugin(null);
	}
	
	/**
	 * Apply this window's workflow to the file associated to this window (through the model) 
	 */
	private void actionApply(){
		
		String path= dialogSave();
		if ("".equals(path)) {
			return;
		}
		
		// Reopen the file only if the data was resized
		TomoData originalModel=getModel();
		
		if(getModel().isResized()){
			originalModel=new TomoData(getModel().getFilePath());
			originalModel.addPropertyChangeListener(this);
			try{
				setReloadingFile(true);
				(new Thread(new ImportDataThread(originalModel,false))).start();
				originalModel.waitForLastImage();
			}catch (Exception ex){
				Xmipp_Tomo.debug("actionApply - unexpected exception", ex);
			}
			
			// iterate through the user actions that make sense
			for (UserAction currentAction : Xmipp_Tomo.getWorkflow(getLastAction())){
				if(currentAction.isNeededForFile()){
					// Xmipp_Tomo.debug("Applying " + currentAction.toString());
					setStatus("Applying " + currentAction.getCommand());
					currentAction.getPlugin().run(originalModel.getImage());
				}
			}
		}
		
		// write 
		saveFile(originalModel, path);
	}
	
	private void captureDialog(){
		Window [] windows=Window.getWindows();
		for(Window w : windows){
			String className= w.getClass().toString();
			// Xmipp_Tomo.debug(className);
			if(className.equals("class ij.gui.GenericDialog")){
				((GenericDialog) w).addDialogListener(this);
				/*for(Component c : w.getComponents()){
					String componentClass = c.getClass().toString();
					Xmipp_Tomo.debug(c.getName());
//					
				}*/
			}
		}
	}
	
	/**
	 * @see DialogListener.dialogItemChanged(GenericDialog gd, AWTEvent e)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e){
		//Xmipp_Tomo.debug("TomoWindow.dialogItemChanged");
		if(gd != null){
			if(getPlugin() != null)
				getPlugin().collectParameters(gd);
		}
		// true so notify continues
		return true;
	}
	
	/* General window events ----------------------------------------------- */
	public void windowClosed(WindowEvent e) {}
	public void windowDeactivated(WindowEvent e) {}
	public void focusLost(FocusEvent e) {}
	public void windowDeiconified(WindowEvent e) {}
	public void windowIconified(WindowEvent e) {}	
	public void windowOpened(WindowEvent e) {}
	public void windowActivated(WindowEvent e) {
	}

	/**
	 * housekeepin'
	 * @param e windowclosing event (unused)
	 */	
	public void windowClosing(WindowEvent e) {
		if (closed)
			return;
		if(isChangeSaved()){
			Xmipp_Tomo.ExitValues choice= dialogYesNoCancel("Close window", "Are you sure you want to close this window?");
			switch (choice){
				case NO:
				case CANCEL:
					return;
			}
		}else{
			Xmipp_Tomo.ExitValues choice= dialogYesNoCancel("Save changes", "Do you want to save changes?");
			switch (choice){
				case YES:
						String path= dialogSave();
						if("".equals(path))
							return;
						saveFile(getModel(),path);
				break;
				case CANCEL:
					return;
			}
		}
		
		// close the window (no return back...)
		setVisible(false);
		dispose();
		WindowManager.removeWindow(this);
		closed=true;
	}

	/* Mouse events ----------------------------------------------- */
	
	// based on ij.gui.ImageCanvas
	public void mouseMoved(MouseEvent e){
		// update status when moving more that 12 pixels
		// if(getCursorDistance(e.getX(),e.getY())> 144)
			updateStatusText();
		setCursorLocation(e.getX(),e.getY());
		
	}
	
	public void mouseExited(MouseEvent e){}
	public void mouseDragged(MouseEvent e){}
	public void mouseClicked(MouseEvent e) {}
	public void mouseEntered(MouseEvent e) {}
	public void mouseReleased(MouseEvent e) {}
	public void mousePressed(MouseEvent e) {}

	// Property changes
	
	public void propertyChange(PropertyChangeEvent event){
		if(TomoData.Properties.NUMBER_OF_PROJECTIONS.name().equals(event.getPropertyName()))
			if(isReloadingFile())
				setStatus("Recovering original file... " + getStatusChar((Integer) (event.getNewValue())));
			else
				setNumberOfProjections((Integer) (event.getNewValue()) ); 
	}
	
	/* Getters/ setters  --------------------------------------------------- */
	
	/**
	 * @return the model
	 */
	public TomoData getModel() {
		return model;
	}

	/**
	 * @param model the model to set
	 */
	public void setModel(TomoData model) {
		this.model = model;
		
		if(this.model != null)
 			// model change events
 			this.model.addPropertyChangeListener(this);		

	}
	
	private JLabel getStatusLabel(){
		return statusLabel;
	}
	
	private char getStatusChar(int i){
		String animation = "\\|/-\\|/-";
		return animation.charAt(i % animation.length());
	}
	
	void setCursorLocation(int x,int y){
		cursorLocation.setLocation(x,y);
	}
	
	
	/**
	 * Update all window components related to the number of projections. The number itself is
	 * stored in the model
	 * @param n - see Tomodata.numberOfProjections
	 */
	public void setNumberOfProjections(int n){
		if(projectionScrollbar != null)
			projectionScrollbar.setMaximum(n);
	}
	
	private int getCursorX(){
		return cursorLocation.x;
	}

	private int getCursorY(){
		return cursorLocation.y;
	}
	
	private double getCursorValue(){
		return getModel().getPixelValue(getCursorX(),getCursorY());
	}
	
	/*
	 * It requires that the model has been set (with setModel() )
	 */
	public String getCursorValueAsString(){
		if(getModel() == null)
			return "";
		
		String cursorValue = String.valueOf(getCursorValue());
		// include the decimal point in the count
		int limit = Math.min(cursorValue.length(), MAX_DIGITS+1);
		return cursorValue.substring(0, limit);
	}

	public String getTitle(){
		String title= TITLE;
		if(getModel() != null)
			title = title + " > " + getModel().getFileName() + " (" + getModel().getWidth() + "x" +  getModel().getHeight() + ")";
		return title; 
	}
	
	/** Show a file open browser and... 
	 * @return the path of the file chosen by the user, or empty string (not null) if the user cancels the dialog
	 */
	private String dialogOpen(){
		OpenDialog od = new OpenDialog("Import file",null);
        String directory = od.getDirectory();
		String fileName = od.getFileName();
		String path= "";
		if(fileName != null)
			path=directory + fileName;
		return path;
	}
	
	/** Show a file save browser and... 
	 * @return the path of the file chosen by the user, or empty string (not null) if the user cancels the dialog
	 */
	private String dialogSave(){
		SaveDialog sd = new SaveDialog("Save as...", getModel().getFilePath(), ".mrc");
	    String directory = sd.getDirectory();
		String fileName = sd.getFileName();
		String path= "";
		if(fileName != null)
			path=directory + fileName;
		return path;
	}
	
	/**
	 * @param title Dialog title
	 * @param message Dialog message
	 * @return Xmipp_Tomo.ExitValues.CANCEL / YES / NO
	 */
	private Xmipp_Tomo.ExitValues dialogYesNoCancel(String title, String message){
		GenericDialog gd = new GenericDialog(title);

		gd.addMessage(message);
        gd.enableYesNoCancel();
        gd.showDialog();
        if (gd.wasCanceled())
           return Xmipp_Tomo.ExitValues.CANCEL;
        else if (gd.wasOKed())
        	return Xmipp_Tomo.ExitValues.YES;
        else
        	return Xmipp_Tomo.ExitValues.NO;
	}
	
	private void saveFile(TomoData model,String path){
		setStatus("Saving...");
		model.setFile(path);
		new TiltSeriesOpener(false).write(model);
		setStatus("Done");
		setChangeSaved(true);
	}
	
    /********************** Main - class testing **************************/

	private void gridBagLayoutTest(){
        Container pane=getContentPane();
        
		JButton button;
    	pane.setLayout(new GridBagLayout());
    	GridBagConstraints c = new GridBagConstraints();
    	
    	button = new JButton("Long-Named Button 4");
    	c.fill = GridBagConstraints.HORIZONTAL;
    	c.ipady = 40;      //make this component tall
    	c.weightx = 0.0;
    	c.gridwidth = 3;
    	c.gridx = 0;
    	c.gridy = 0;
    	pane.add(button, c);

    	setVisible(true);
	}
	
	private void layoutTest(){
		addMainPanels();
		setVisible(true);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
	    TomoWindow w= new TomoWindow();
	    /* frame.setSize(400, 400);
	    frame.setLocationRelativeTo(null);
	    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); */

	    
	    //w.gridBagLayoutTest();
	    w.layoutTest();
	    
	}

	/**
	 * @return the firstAction
	 */
	public UserAction getFirstAction() {
		return firstAction;
	}

	/**
	 * @param firstAction the firstAction to set
	 */
	public void setFirstAction(UserAction firstAction) {
		this.firstAction = firstAction;
	}

	/**
	 * @return the lastAction
	 */
	public UserAction getLastAction() {
		if(lastAction==null)
			lastAction=firstAction;
		return lastAction;
	}

	/**
	 * @param lastAction the lastAction to set
	 */
	public void setLastAction(UserAction lastAction) {
		this.lastAction = lastAction;
	}

	/**
	 * @return the windowId
	 */
	public int getWindowId() {
		return windowId;
	}

	/**
	 * @param windowId the windowId to set
	 */
	public void setWindowId(int windowId) {
		this.windowId = windowId;
	}

	/**
	 * @return the plugin
	 */
	public Plugin getPlugin() {
		return plugin;
	}

	/**
	 * @param plugin the plugin to set
	 */
	public void setPlugin(Plugin plugin) {
		this.plugin = plugin;
	}

	/**
	 * @return the reloadingFile
	 */
	private boolean isReloadingFile() {
		return reloadingFile;
	}

	/**
	 * @param reloadingFile the reloadingFile to set
	 */
	private void setReloadingFile(boolean reloadingFile) {
		this.reloadingFile = reloadingFile;
	}

	/**
	 * @return the changeSaved
	 */
	private boolean isChangeSaved() {
		return changeSaved;
	}

	/**
	 * @param changeSaved the changeSaved to set
	 */
	private void setChangeSaved(boolean changeSaved) {
		this.changeSaved = changeSaved;
	}
}
