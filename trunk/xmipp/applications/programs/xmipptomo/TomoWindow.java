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

import ij.WindowManager;
import ij.gui.ImageCanvas;
import ij.gui.ImageLayout;

import java.awt.*;
import java.awt.event.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JScrollBar;
import javax.swing.JTextField;



/**
 * @author jcuenca
 * Custom visualization window, using Swing and MVC paradign
 */
public class TomoWindow extends JFrame implements WindowListener, AdjustmentListener,MouseMotionListener, PropertyChangeListener{
	// for serialization only
	private static final long serialVersionUID = -4063711975454855701L;
	
	// total number of  digits to display as a String
	private static int MAX_DIGITS = 9;
	// for scrollbars
	private static int extent=20;

	// Window title (without extra data; see getTitle() )
	final static String TITLE="XmippTomo";
	private boolean closed=true;
	
	/* Window components */
	/* 3 main panels */
	Container viewPanel,controlPanel,statusPanel;
	// Text fields and labels
	JTextField tiltTextField;
	JLabel statusLabel;
	
	// The canvas where current projection is displayed
	ImageCanvas ic;
	// Control scrollbars
	// JScrollBar projectionScrollbar;
	LabelScrollbar projectionScrollbar;
	
	// the model stores the data that this window shows
	private TomoData model;
	
	Point cursorLocation=new Point();
	
	/* CONSTRUCTORS -------------------------------------------------------- */
	
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

	/* METHODS ------------------------------------------------------------- */
	
	/**
	 * Create all window components and associate each text field with model's documents
	 * @param model
	 */
	public void create(TomoData model){
				
		setModel(model);
		
		// set this window as listener of general keyboard&mouse events and model change events
 		addWindowListener(this);
 		model.addPropertyChangeListener(this);
 		
 		// EXIT_ON_CLOSE finishes ImageJ too...
		
 		// Window title
 		setTitle(getTitle());
 		
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
		ic.addMouseMotionListener(this);
		
		// CONTROLS PANEL
		
		controlPanel = new Container();
		FlowLayout fl2=new FlowLayout();
		fl2.setAlignment(FlowLayout.CENTER);
		controlPanel.setLayout(fl2);
		constraints.gridx=0; constraints.gridy = 1;
		constraints.fill = GridBagConstraints.HORIZONTAL;		
		gb1.setConstraints(controlPanel, constraints);
		getContentPane().add(controlPanel);
		JLabel tiltTextLabel=new JLabel("Tilt");
		tiltTextField=new JTextField(5);
		controlPanel.add(tiltTextLabel);		
		controlPanel.add(tiltTextField);
		getModel().setTiltModel(tiltTextField.getDocument());
		
		projectionScrollbar= new LabelScrollbar(1,getModel().getNumberOfProjections());
		projectionScrollbar.setText("Projection #");
		projectionScrollbar.addAdjustmentListener(this);
		controlPanel.add(projectionScrollbar);
		
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
		updateStatusText();
		
		closed=false;
	}
	
	/**
	 * warning: naming this method "show" leads to infinite recursion (since setVisible calls show in turn)
	 */
	public void display(){
		 pack();
	     setVisible(true);
	}
	
	void updateStatusText(){
		// current/total projections are shown in the projection scrollbar itself
		getStatusLabel().setText("x = " + getCursorX() + ", y = " + getCursorY() + ", value = " + getCursorValueAsString());
	}
	
	int getCursorDistance(int x, int y){
		return (int) cursorLocation.distance(x,y);
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
			ic.setImageUpdated();
			ic.repaint();
			updateStatusText(); // projection number
		}
		notify();
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
	
	int getCursorX(){
		return cursorLocation.x;
	}

	int getCursorY(){
		return cursorLocation.y;
	}
	
	double getCursorValue(){
		return getModel().getPixelValue(getCursorX(),getCursorY());
	}
	
	public String getCursorValueAsString(){
		String cursorValue = String.valueOf(getCursorValue());
		// include the decimal point in the count
		int limit = Math.min(cursorValue.length(), MAX_DIGITS+1);
		return cursorValue.substring(0, limit);
	}

	public String getTitle(){
		return TITLE + " > " + getModel().getFileName() + " (" + getModel().getWidth() + "x" +  getModel().getHeight() + ")";
	}
}
