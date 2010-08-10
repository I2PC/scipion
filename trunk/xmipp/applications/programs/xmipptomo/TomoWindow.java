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
import ij.gui.ImageCanvas;
import ij.gui.ImageLayout;
import ij.io.OpenDialog;

import java.awt.*;
import java.awt.event.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
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
 * Custom visualization window, using Swing and MVC paradign
 */
/** Protocol for adding new buttons to the workflow/UI:
 *  - set the label in Commands enum as static String,
 *  - run addButton in Xmipp_Tomo constructor, 
 *  - update actionPerformed method
 */
public class TomoWindow extends JFrame implements WindowListener, AdjustmentListener,MouseMotionListener,ActionListener, PropertyChangeListener{
	// for serialization only
	private static final long serialVersionUID = -4063711975454855701L;
	
	// total number of  digits to display as a String
	private static int MAX_DIGITS = 9;
	// Window title (without extra data; see getTitle() )
	private final static String TITLE="XmippTomo";
	private static int MENUPANEL_MINWIDTH=400, MENUPANEL_MINHEIGHT=100;
	
	// labels of buttons
	private static enum ButtonLabels {
		CMD_LOAD("Load",""),CMD_SAVE("Save",""),
		CMD_GAUSSIAN("Gaussian","Gaussian Blur..."),CMD_MEDIAN("Median", "Median..."),CMD_SUB_BACKGROUND("Bandpass - substract background","Bandpass Filter...");

		private final String label;
		private final String imageJCmd;
		
		ButtonLabels(String l,String cmd) { this.label=l;this.imageJCmd=cmd;}
		public String label(){return label;}
		public String imageJCmd(){return imageJCmd;}
	};

	private boolean closed=true;
	
	/* Window components */
	/* 4 main panels */
	private JTabbedPane menuPanel;
	private JPanel viewPanel,controlPanel,statusPanel;
	
	// Text fields and labels
	private JTextField tiltTextField;
	private JLabel statusLabel;
	
	// The canvas where current projection is displayed
	private ImageCanvas ic;
	
	// Control scrollbars
	private LabelScrollbar projectionScrollbar;
	
	// the model stores the data that this window shows
	private TomoData model=null;
	
	Point cursorLocation=new Point();
	
	// implement UI concurrency with private class that extends Thread
	// ImportDataThread: load & display projections in parallel
	private class ImportDataThread extends Thread{
		private TomoData dataModel;
		private String dataPath;
		
		ImportDataThread(TomoData model,String path){
			dataModel=model;
			dataPath=path;
		}
		
		public void run() {
			try{
				dataModel.import_data(dataPath);
			}catch (IOException ex){
				Xmipp_Tomo.debug("ImportDataThread.run - Error opening file");
			}catch (InterruptedException ex){
				Xmipp_Tomo.debug("ImportDataThread.run - Interrupted exception");
			}catch (Exception ex){
				Xmipp_Tomo.debug("ImportDataThread.run - unexpected exception", ex);
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
		WindowManager.addWindow(this);
		addMainPanels();
 		// EXIT_ON_CLOSE finishes ImageJ too...
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

		// BoxLayout component order is sequential, so all the panels must be added in order (even if they are empty)
 		
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
		pack();
	}
	
	/**
	 * It requires that the model has been set (with setModel() ), and that at least 1 image has been loaded
	 * @param model
	 */
	private void addView(){
		if(getModel() == null)
			return;
		
		ic=new ImageCanvas(getModel().getImage());
		/* iw=new ImageWindow(getModel().getImage());
		WindowManager.setCurrentWindow(iw);
		iw.setVisible(false);
		iw.hide();
		Xmipp_Tomo.debug("Visible: " + iw.isVisible());
*/
		viewPanel.setLayout(new ImageLayout(ic));
		viewPanel.add(ic);
		ic.addMouseMotionListener(this);
	
		pack();// adjust window size to host the new view panel
	}
	
	/**
	 * It also requires that the model has been set (with setModel() ), and that basic info (number of projections,...) has been loaded
	 * @param model
	 */
	private void addControls(){
		if(getModel() == null)
			return;
			

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
		addButton(ButtonLabels.CMD_LOAD.label(),fileTabPanel);
		addButton(ButtonLabels.CMD_SAVE.label(),fileTabPanel);
		menuPanel.addTab("File",fileTabPanel); // try to add the tab when it's ready (all controls inside)

		// Preprocessing
		JPanel preprocTabPanel=new JPanel();
		addButton(ButtonLabels.CMD_GAUSSIAN.label(),preprocTabPanel);
		addButton(ButtonLabels.CMD_MEDIAN.label(),preprocTabPanel);
		addButton(ButtonLabels.CMD_SUB_BACKGROUND.label(),preprocTabPanel);
		menuPanel.addTab("Preprocessing",preprocTabPanel);
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
		ic.setImageUpdated();
		ic.repaint();
	}
	
	/**
	 * Scrollbar events
	 * @param e
	 */
	public synchronized void adjustmentValueChanged(AdjustmentEvent e){
		// Xmipp_Tomo.debug("TomoWindow.adjustmentValueChanged"+e.getSource().toString());
		if(e.getSource()==projectionScrollbar){
			// scrollbar range is 1..N, while projections array is 0..N-1
			getModel().setCurrentProjection(projectionScrollbar.getValue()-1);
			refreshImageCanvas();
			updateStatusText(); // projection number
		}
		notify();
	}
	
	/* button/keyboard events
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	public void actionPerformed(ActionEvent e) {
	
		String label = e.getActionCommand();

		// select proper action method based on button's label
		if (label==null)
			return;
		else if (label.equals(ButtonLabels.CMD_LOAD.label())){
			this.actionLoad();
		}else if (label.equals(ButtonLabels.CMD_GAUSSIAN.label())){
			this.actionRunIjCmd(ButtonLabels.CMD_GAUSSIAN.label(),ButtonLabels.CMD_GAUSSIAN.imageJCmd());
		}else if (label.equals(ButtonLabels.CMD_MEDIAN.label())){
			this.actionRunIjCmd(ButtonLabels.CMD_MEDIAN.label(),ButtonLabels.CMD_MEDIAN.imageJCmd());
		}else if (label.equals(ButtonLabels.CMD_SUB_BACKGROUND.label())){
			this.actionRunIjCmd(ButtonLabels.CMD_SUB_BACKGROUND.label(),ButtonLabels.CMD_SUB_BACKGROUND.imageJCmd());
		}
	} // actionPerformed end
	
	private void actionLoad(){
		String path = browseFile();
		
		try{
			// import data in one thread and hold this thread until the first projection is loaded
			if(getModel() == null)
				setModel(new TomoData());
			
			(new Thread(new ImportDataThread(getModel(),path))).start();
			getModel().waitForFirstImage();
		}catch (InterruptedException ex){
			Xmipp_Tomo.debug("Xmipp_Tomo - Interrupted exception");
		}

		addView();
		addControls();
	}
	
	private void actionRunIjCmd(String label,String cmd){
		// set current IJ window (windowmanager)
		if(label != null)
			setStatus(label + " - started");
		
		WindowManager.setTempCurrentImage(getModel().getImage());
		// IJ.doCommand(ButtonLabels.CMD_GAUSSIAN.imageJCmd());
		IJ.run(cmd);
		
		refreshImageCanvas();
		
		setStatus("Done");
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

		setVisible(false);
		dispose();
		WindowManager.removeWindow(this);
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
	
	/** Show a file browser and... 
	 * @return the path of the file chosen by the user
	 */
	private String browseFile(){
		OpenDialog od = new OpenDialog("Import file",null);
        String directory = od.getDirectory();
		String fileName = od.getFileName();
		String path= directory + fileName;
		return path;
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
}
