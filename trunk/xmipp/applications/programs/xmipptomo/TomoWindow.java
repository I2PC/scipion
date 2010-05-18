package xmipptomo;

import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ImageCanvas;
import ij.gui.ImageLayout;


import java.awt.*;
import java.awt.event.*;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JScrollBar;
import javax.swing.JTextField;


/**
 * Custom X-Y-Z visualization window, using Swing
 */

/**
 * @author jcuenca
 *
 */
public class TomoWindow extends JFrame implements WindowListener, AdjustmentListener{
	private static final long serialVersionUID = -4063711975454855701L;

	final static String TITLE="XmippTomo XYZ";
	private boolean closed=true;
	
	/* Window components */
	/* 3 main panels */
	Container viewPanel,controlPanel,statusPanel;
	JTextField tiltTextField;
	JLabel statusLabel;
	ImageCanvas ic;
	JScrollBar xScrollbar, yScrollbar, zScrollbar, tiltScrollbar;
	
	private TomoData model;
	
	/**
	 * @throws HeadlessException
	 */
	public TomoWindow() throws HeadlessException {
		// TODO Auto-generated constructor stub
		super(TITLE);
	}

	/**
	 * @param gc
	 */
	public TomoWindow(GraphicsConfiguration gc) {
		super(gc);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param title
	 * @throws HeadlessException
	 */
	public TomoWindow(String title) throws HeadlessException {
		super(title);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param title
	 * @param gc
	 */
	public TomoWindow(String title, GraphicsConfiguration gc) {
		super(title, gc);
		// TODO Auto-generated constructor stub
	}
	
	/**
	 * Associate each text field with model's documents
	 * @param imp
	 * @param ti
	 */
	public void create(TomoData model){
		
		JLabel tiltTextLabel;
		
		setModel(model);
		
 		addWindowListener(this);
 		// EXIT_ON_CLOSE finishes ImageJ too...
 		// setDefaultCloseOperation(EXIT_ON_CLOSE);
 		// getContentPane().add(new JLabel("Hello"));
 		
 	
		GridBagLayout gb1= new GridBagLayout();
		getContentPane().setLayout(gb1);
		GridBagConstraints constraints=new GridBagConstraints();
		
		// VIEW PANEL
		// View panel should fill all available X & Y space, the others should grow 
		// only horizontally
		constraints.fill = GridBagConstraints.BOTH;

		viewPanel = new Container();
		constraints.gridx=0; constraints.gridy = 0;
		gb1.setConstraints(viewPanel, constraints);
		getContentPane().add(viewPanel);
		
		
		ic=new ImageCanvas(getModel().getImage());
		viewPanel.setLayout(new ImageLayout(ic));
		viewPanel.add(ic);
		
		
		// CONTROLS PANEL
		controlPanel = new Container();
		FlowLayout fl2=new FlowLayout();
		fl2.setAlignment(FlowLayout.CENTER);
		controlPanel.setLayout(fl2);
		constraints.gridx=0; constraints.gridy = 1;
		constraints.fill = GridBagConstraints.HORIZONTAL;		
		gb1.setConstraints(controlPanel, constraints);
		getContentPane().add(controlPanel);
		tiltTextLabel=new JLabel("Tilt");
		tiltTextField=new JTextField(5);
		controlPanel.add(tiltTextLabel);		
		controlPanel.add(tiltTextField);
		getModel().setTiltModel(tiltTextField.getDocument());
		//tiltTextField.setText(Integer.toString(0));
		// X Y Z scrollbars
		tiltScrollbar= new JScrollBar(JScrollBar.HORIZONTAL);
		// Calculate tilt range - not really needed to select slice
		// float min=TomoInfo.getMinTilt(ti), max=TomoInfo.getMaxTilt(ti);
		tiltScrollbar.setMinimum(0);
		tiltScrollbar.setMaximum(getModel().getNSlices());
		tiltScrollbar.addAdjustmentListener(this);
		controlPanel.add(tiltScrollbar);
		
		
		
		// STATUS PANEL
		statusPanel = new Container();
		FlowLayout fl=new FlowLayout();
		fl.setAlignment(FlowLayout.LEFT);
		statusPanel.setLayout(fl);
		constraints.gridx=0; constraints.gridy = 2;
		constraints.fill = GridBagConstraints.HORIZONTAL;
		gb1.setConstraints(statusPanel, constraints);
		getContentPane().add(statusPanel);
		statusLabel=new JLabel("Done");
		statusPanel.add(statusLabel);
		
		closed=false;
	
	}
	
	/**
	 * naming this show leads to infinite recursion since setVisible calls show in turn
	 */
	public void display(){
		 pack();
	     setVisible(true);
	}
	
	/**
	 * @param e
	 *  Scrollbar events
	 */
	public synchronized void adjustmentValueChanged(AdjustmentEvent e){
		if(e.getSource()==tiltScrollbar){
			// imp.updatePosition(1, tiltScrollbar.getValue() +1, 1);
			getModel().setCurrentSlice(tiltScrollbar.getValue()+1);
			ic.setImageUpdated();
			ic.repaint();
			// imp.updateAndDraw();
			// int tiltValue = tiltScrollbar.getValue();
			// int tiltValue = imp.getSlice();
			//float tiltValue = 
			// tiltTextField.setText(Integer.toString(tiltValue));
		}
		notify();
	}
	
	/* General events */
	public void windowClosed(WindowEvent e) {}
	public void windowDeactivated(WindowEvent e) {}
	public void focusLost(FocusEvent e) {}
	public void windowDeiconified(WindowEvent e) {}
	public void windowIconified(WindowEvent e) {}	
	public void windowOpened(WindowEvent e) {}

	public void windowActivated(WindowEvent e) {
		// WindowManager.setCurrentWindow(this);
	}
	
	public void windowClosing(WindowEvent e) {
		if (closed)
			return;

/*		WindowManager.setCurrentWindow(this);
			IJ.doCommand("Close"); */
 			setVisible(false);
			dispose();
			WindowManager.removeWindow(this);
	}

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
	}

}
