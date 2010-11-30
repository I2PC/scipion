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

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.ImageLayout;
import ij.gui.ImageWindow;
import ij.io.SaveDialog;

import java.awt.*;
import java.awt.event.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Hashtable;
import java.util.LinkedList;

import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.Timer;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.plaf.FontUIResource;

/**
 * @author jcuenca
 * Custom visualization window, using Swing and MVC paradigm
 * 
 * implements: WindowListener, AdjustmentListener (scrollbars), MouseMotionListener, ActionListener, PropertyChangeListener to handle GUI events
 * 			DialogListener to capture GenericDialogs, ComponentListener to handle window resizing
 */
/**
 * Protocol for adding new buttons to the workflow/UI:
 * - declare its command in Command
 * - define the action method in TomoController
 * - include the command in its menu list (for example, commandsMenuFile) 
 * - set a proper enable/disable cycle along user interaction
 * (tipically start disabled, enable it when user can click it)
 */
public class TomoWindow extends ImageWindow implements WindowListener,
		AdjustmentListener, MouseMotionListener, ActionListener,
		PropertyChangeListener, DialogListener, ComponentListener {
	// for serialization only - just in case
	private static final long serialVersionUID = -4063711975454855701L;

	// total number of digits to display as a String
	private static int MAX_DIGITS = 9;
	// Window title (without extra data; see getTitle() )
	private final static String TITLE = "XmippTomo";
	// minimum dimension of tabbed menu panel (window top)
	private static int MENUPANEL_MINWIDTH = 400, MENUPANEL_MINHEIGHT = 100;
	// maximum window size
	private static Rectangle maxWindowSize = new Rectangle(800, 800);
	// miliseconds to wait for a Plugin dialog to display (so the dialog capture
	// can start)
	public static int WAIT_FOR_DIALOG = 800;

	public static String PLAY_ICON = "resources/icon-play.png",
			PAUSE_ICON = "resources/icon-pause.png";
	
	private Hashtable<String, ImageIcon> icons = new Hashtable<String, ImageIcon>();

	private Hashtable<String, JButton> buttons = new Hashtable<String, JButton>();

	// Lists of commands for each menu tab
	private static java.util.List<Command> commandsMenuFile = new LinkedList<Command>() {
		{
			add(Command.LOAD);
			add(Command.XRAY);
			add(Command.SAVE);
			add(Command.DEFINE_TILT);
			add(Command.MEASURE);
		}
	};

	private static java.util.List<Command> commandsMenuPreproc = new LinkedList<Command>() {
		{
			add(Command.GAUSSIAN);
			add(Command.MEDIAN);
			add(Command.SUB_BACKGROUND);
			add(Command.ENHANCE_CONTRAST);
			add(Command.APPLY);
		}
	};

	private static java.util.List<Command> commandsMenuAlign = new LinkedList<Command>() {
		{
			add(Command.ALIGN_AUTO);
			add(Command.ALIGN_MANUAL);
			add(Command.ALIGN_CORRELATION);
		}
	};

	private static java.util.List<Command> commandsMenuDebug = new LinkedList<Command>() {
		{
			add(Command.PRINT_WORKFLOW);
		}
	};

	// which commands to enable after loading a file
	private static java.util.List<Command> commandsEnabledAfterLoad = new LinkedList<Command>() {
		{
			add(Command.SAVE);
			add(Command.DEFINE_TILT);
			add(Command.GAUSSIAN);
			add(Command.MEDIAN);
			add(Command.SUB_BACKGROUND);
			add(Command.ENHANCE_CONTRAST);
			add(Command.APPLY);
			add(Command.PRINT_WORKFLOW);
			add(Command.MEASURE);
			add(Command.PLAY);
		}
	};

	// Menu labels
	private enum Menu {
		FILE("File"), PREPROC("Preprocess"), ALIGN("Align"), RECONSTRUCTION("Reconstruction"), DEBUG("Debug");

		private String label;

		private Menu(String l) {
			this.label = l;
		}

		public final String toString() {
			return label;
		}
	}

	private boolean closed = true;
	private Command.State lastCommandState=Command.State.IDLE;
	// true if changes are saved (so you can close the window without pain)
	private boolean changeSaved = true;
	// true if current projection is changing automatically
	private boolean playing = false;
	// window identifier - useful when you've got more than 1 window
	private int windowId = -1;
	private TomoController controller;

	/* Window GUI components */
	/* 4 main panels */
	private JTabbedPane menuPanel;
	private JPanel viewPanel, controlPanel, statusPanel;
	private JFrame realWindow;

	// Text fields and labels
	private JTextField tiltTextField;
	private JLabel statusLabel;

	// Control scrollbars
	private LabelScrollbar projectionScrollbar;

	// the model stores the data that this window shows - @see TomoData
	private TomoData model = null;

	Point cursorLocation = new Point();

	// action history hints - so we can navigate the workflow
	private UserAction firstAction, lastAction;

	// for window resizing
	Timer timer;

	// current plugin - DialogCapture thread uses it (when dialotItemChanged is
	// called)
	private Plugin plugin;

	/* METHODS ------------------------------------------------------------- */

	/*
	 * For the buggy GridBagLayout... // row: 1..N private void
	 * addToGrid(Container c,int row, int constraints_fill){ GridBagConstraints
	 * constraints=new GridBagConstraints(); constraints.gridx=0;
	 * constraints.gridy = row-1; constraints.anchor = GridBagConstraints.SOUTH;
	 * constraints.fill = constraints_fill; //gbl.setConstraints(c,
	 * constraints); getContentPane().add(c,constraints); }
	 */

	/**
	 * Create a window with its main panels (but no data/model available) and
	 * show it
	 */
	public TomoWindow() {
		super("Hola");

		setController(new TomoController(this));
		
		// change default font
		Font defaultBoldFont = new Font("Dialog", Font.BOLD, 14), defaultFont = new Font(
				"Dialog", Font.PLAIN, 14);
		UIManager.put("Label.font", new FontUIResource(defaultBoldFont));
		UIManager.put("Button.font", new FontUIResource(defaultBoldFont));
		UIManager.put("Panel.font", new FontUIResource(defaultBoldFont));
		UIManager.put("TextField.font", new FontUIResource(defaultFont));
		UIManager.put("TabbedPane.font", new FontUIResource(defaultBoldFont));

		// turn on anti-aliasing for smooth fonts.
		// System.setProperty( "swing.aatext", "true" );
		// ...or (Java 1.6)
		System.setProperty("swing.useSystemAAFontSettings", "lcd");
		
		// disable dynamic resizing - it's OS dependent...
		// Toolkit.getDefaultToolkit().setDynamicLayout(false);
		// set up a fake ImagePlus, while we get to load the real one
		imp = new ImagePlus();
		imp.setWindow(this);
		realWindow = new JFrame(TITLE);
		// set this window as listener of general keyboard&mouse events
		realWindow.addWindowListener(this);
		realWindow.addComponentListener(this);
		// EXIT_ON_CLOSE finishes ImageJ too...
		// DO NOTHING allows for closing confirmation dialogs and the like
		// setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		WindowManager.addWindow(this);
		addMainPanels();
	}

	public TomoWindow(int windowId) {
		this();
		setWindowId(windowId);
	}

	// JFrame methods
	private Container getContentPane() {
		return realWindow.getContentPane();
	}

	private void addButton(Command cmd, Container panel) throws Exception {
		TomoAction a = new TomoAction(getController(), cmd);
		JButton b = new JButton(a);
		panel.add(b);
		buttons.put(cmd.getId(), b);
	}

	public JButton getButton(String id) {
		return buttons.get(id);
	}

	/**
	 * Add the main panels with their components, and associate each text field
	 * with model's documents (if a model is already set)
	 */
	private void addMainPanels() {
		// Window title
		realWindow.setTitle(getTitle());

		realWindow.setLayout(new BoxLayout(getContentPane(),
				BoxLayout.PAGE_AXIS));

		// BoxLayout component order is sequential, so all the panels
		// MUST be added in order (even if they are empty)

		// TABBED MENU PANEL
		menuPanel = new JTabbedPane();
		try{
			addMenuTabs();
		}catch (Exception ex){
			Xmipp_Tomo.debug("Tomowindow.addMainPanels - addmenutabs " + ex.toString());
		}
		menuPanel.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
		menuPanel.setPreferredSize(new Dimension(MENUPANEL_MINWIDTH,
				MENUPANEL_MINHEIGHT));
		getContentPane().add(menuPanel);

		// VIEW & CONTROLS PANELS
		viewPanel = new JPanel();
		getContentPane().add(viewPanel);

		controlPanel = new JPanel();
		getContentPane().add(controlPanel);

		if (getModel() != null) {
			addView();
			addControls();
		}

		// STATUS PANEL
		statusPanel = new JPanel();
		FlowLayout fl = new FlowLayout();
		fl.setAlignment(FlowLayout.LEFT);
		statusPanel.setLayout(fl);
		statusPanel.setMaximumSize(new Dimension(Short.MAX_VALUE, 50));
		statusPanel.setBorder(BorderFactory.createMatteBorder(1, 0, 0, 0, Color.BLACK));

		statusLabel = new JLabel("Ready");
		statusPanel.add(statusLabel);
		updateStatusText();
		getContentPane().add(statusPanel);

		closed = false;
		realWindow.pack(); // adjust GUI size
	}

	/**
	 * It requires that the model has been set (with setModel() ), and that at
	 * least 1 image has been loaded
	 * 
	 * @param model
	 */
	public void addView() {
		if (getModel() == null)
			return;

		// remove previous canvas if any
		for (int i = 0; i < viewPanel.getComponentCount(); i++)
			viewPanel.remove(i);

		setCanvas(new ImageCanvas(getModel().getImage()));
		viewPanel.setLayout(new ImageLayout(getCanvas()));
		viewPanel.add(getCanvas());
		getCanvas().addMouseMotionListener(this);
		addPropertyChangeListener(this);

		realWindow.pack(); // adjust window size to host the new view panel
	}

	/**
	 * It also requires that the model has been set (with setModel() ), and that
	 * basic info (number of projections,...) has been loaded
	 * 
	 * @param model
	 */
	public void addControls() {
		// only add controls if they are not available already
		if (projectionScrollbar != null){
			projectionScrollbar.setMaximum(getModel().getNumberOfProjections());
			getModel().setTiltModel(tiltTextField.getDocument());
			return;
		}
		
		/* remove previous controls if any
		for (int i = 0; i < controlPanel.getComponentCount(); i++)
			controlPanel.remove(i); */

		FlowLayout fl2 = new FlowLayout();
		fl2.setAlignment(FlowLayout.CENTER);
		controlPanel.setLayout(fl2);

		JLabel tiltTextLabel = new JLabel("Tilt");
		tiltTextField = new JTextField(5);
		// the user won't change angles one by one - he will change all of them
		// at once with the tilt change dialog
		tiltTextField.setEditable(false);
		controlPanel.add(tiltTextLabel);
		controlPanel.add(tiltTextField);
		getModel().setTiltModel(tiltTextField.getDocument());

		projectionScrollbar = new LabelScrollbar(1, getModel()
				.getNumberOfProjections());
		projectionScrollbar.setText("Projection #");
		projectionScrollbar.addAdjustmentListener(this);
		controlPanel.add(projectionScrollbar);

		// addButton(PLAY, controlPanel);
		try {
			addButton(Command.PLAY, controlPanel);
		} catch (Exception ex) {
			Xmipp_Tomo.debug("addControls - play");
		}
		realWindow.pack(); // adjust window size to host the new controls panel

	}

	/**
	 * Add menu tabs populated with their buttons (specified in the lists at the
	 * beginning of the class)
	 */
	private void addMenuTabs() throws Exception {
		// File
		JPanel fileTabPanel = new JPanel(false); // false = no double buffer
		// (flicker)
		// menuPanel.setMnemonicAt(0, KeyEvent.VK_1);
		for (Command c : commandsMenuFile)
			addButton(c, fileTabPanel);

		// try to add the tab when it's ready (all controls inside)
		menuPanel.addTab(Menu.FILE.toString(), fileTabPanel);

		// Preprocessing
		JPanel preprocTabPanel = new JPanel();
		for (Command c : commandsMenuPreproc)
			addButton(c, preprocTabPanel);

		menuPanel.addTab(Menu.PREPROC.toString(), preprocTabPanel);

		// Align
		JPanel alignTabPanel = new JPanel();
		for (Command c : commandsMenuAlign)
			addButton(c, alignTabPanel);
		menuPanel.addTab(Menu.ALIGN.toString(), alignTabPanel);
		
		// Reconstruction
		JPanel reconstructionTabPanel = new JPanel();
		/*for (Command c : commandsMenuAlign)
			addButton(c, alignTabPanel);*/
		menuPanel.addTab(Menu.RECONSTRUCTION.toString(), reconstructionTabPanel);
		menuPanel.setEnabledAt(3, false);
		
		// Debugging
		if (Xmipp_Tomo.TESTING == 1) {
			JPanel debugTabPanel = new JPanel();
			for (Command c : commandsMenuDebug)
				addButton(c, debugTabPanel);
			// addButton(PRINT_WORKFLOW,debugTabPanel,false);
			menuPanel.addTab(Menu.DEBUG.toString(), debugTabPanel);
		}
	}

	/**
	 * warning: naming this method "show" leads to infinite recursion (since
	 * setVisible calls show in turn)
	 */
	public void display() {
		// pack();
		setVisible(false);
		realWindow.pack();
		realWindow.setVisible(true);
	}

	public void setStatus(String text) {
		getStatusLabel().setText(text);
	}

	private void updateStatusText() {
		// current/total projections are shown in the projection scrollbar
		// itself
		if (getModel() != null)
			setStatus("x = " + getCursorX() + ", y = " + getCursorY()
					+ ", value = " + getCursorValueAsString());
	}

	int getCursorDistance(int x, int y) {
		return (int) cursorLocation.distance(x, y);
	}

	public void refreshImageCanvas() {
		if (getCanvas() != null) {
			getCanvas().setImageUpdated();
			getCanvas().repaint();
		}
	}

	/**
	 * Scrollbar events
	 * 
	 * @param e
	 */
	public synchronized void adjustmentValueChanged(AdjustmentEvent e) {
		// Xmipp_Tomo.debug("TomoWindow.adjustmentValueChanged"+e.getSource().toString());
		if (e.getSource() == projectionScrollbar) {
			// scrollbar range is 1..N, like projections array
			getModel().setCurrentProjection(projectionScrollbar.getValue());
			refreshImageCanvas();
			updateStatusText(); // projection number
		}
		notify();
	}

	public void enableButtonsAfterLoad() {
		for (Command c : commandsEnabledAfterLoad) {
			getButton(c.getId()).setEnabled(true);
		}
	}

	public void setImagePlusWindow() {
		imp = getModel().getImage();
		imp.setWindow(this);
	}

	/*
	 * button/keyboard events - run in a different thread to not to block the
	 * GUI
	 * 
	 * @see
	 * java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	public void actionPerformed(ActionEvent e) {

		String label = e.getActionCommand();
		if (label == null) {
			// Timer event
			realWindow.pack();
			// Xmipp_Tomo.debug("Timer");
		}
	} // actionPerformed end

	/* -------------------------- Action management -------------------- */

	public void addUserAction(UserAction a) {
		Xmipp_Tomo.addUserAction(getLastAction(), a);
		setLastAction(a);
	}

	public void changeIcon(String buttonId, String iconName) {
		getButton(buttonId).setIcon(getIcon(iconName));
	}

	public void captureDialog() {
		Window[] windows = Window.getWindows();
		for (Window w : windows) {
			String className = w.getClass().toString();
			// Xmipp_Tomo.debug(className);
			if (className.equals("class ij.gui.GenericDialog")) {
				((GenericDialog) w).addDialogListener(this);
				/*
				 * for(Component c : w.getComponents()){ String componentClass =
				 * c.getClass().toString(); Xmipp_Tomo.debug(c.getName()); // }
				 */
			}
		}
	}

	/**
	 * @see DialogListener.dialogItemChanged(GenericDialog gd, AWTEvent e)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		// Xmipp_Tomo.debug("TomoWindow.dialogItemChanged");
		if (gd != null) {
			if (getPlugin() != null)
				getPlugin().collectParameters(gd);
		}
		// true so notify continues
		return true;
	}

	/* Component events ---------------------------------------------------- */
	public void componentHidden(ComponentEvent e) {
	}

	public void componentMoved(ComponentEvent e) {
	}

	// This event is called many times as the user resizes the window - it may
	// be better to capture
	// a propertychange event like "windowresized" (meaning the user released
	// the mouse and finished
	// resizing. Unluckily, mouse release is not fired while resizing a window)
	// if it exists
	public void componentResized(ComponentEvent e) {
		Component c = e.getComponent();
		// Xmipp_Tomo.debug("componentResized event from " +
		// c.getClass().getName());
		// + "; new size: " + c.getSize().width + ", "
		// + c.getSize().height + " - view size:" + viewPanel.getSize());
		getTimer().stop();
		// getTimer().setDelay(2000);
		timer.start();
		resizeView();
	}

	/**
	 * Fit the image canvas to the space available in the view, keeping the
	 * aspect ratio
	 */
	public void resizeView() {
		if ((viewPanel == null) || (getCanvas() == null))
			return;

		double factorWidth = viewPanel.getWidth()
				/ getCanvas().getPreferredSize().getWidth();
		double factorHeight = viewPanel.getHeight()
				/ getCanvas().getPreferredSize().getHeight();
		double factor = 1;

		if (factorWidth < 1) {
			if (factorHeight < 1)
				// apply maximum of both (since both are < 1 )
				factor = Math.max(factorWidth, factorHeight);
			else
				factor = factorWidth;
		} else {
			if (factorHeight < 1)
				factor = factorHeight;
			else
				factor = Math.min(factorWidth, factorHeight);
		}

		// if size change is minimal or none, don't resize canvas
		if (Math.abs(factor - 1) < 0.03)
			return;

		// Xmipp_Tomo.debug("resize: " + factorWidth + "," + factorHeight + ", "
		// + factor );
		resizeView(factor);

	}

	public void resizeView(double factor) {
		double w = getCanvas().getWidth() * factor;
		double h = getCanvas().getHeight() * factor;

		// getCanvas().resizeCanvas((int)w, (int)h);
		/*
		 * ImageCanvas canvas = getModel().getImage().getCanvas();
		 * getCanvas().setMagnification(factor); canvas.setSourceRect(new
		 * Rectangle(0, 0, (int)(w/factor), (int)(h/factor)));
		 * canvas.setDrawingSize((int)w, (int)h); realWindow.pack();
		 * canvas.repaint();
		 */
		WindowManager.setTempCurrentImage(getModel().getImage());
		// Xmipp_Tomo.debug("Set... " + "zoom="+ ((int) (factor * 100)));
		IJ.run("Set... ", "zoom=" + ((int) (factor * 100)));
		refreshImageCanvas();

	}

	public void componentShown(ComponentEvent e) {
	}

	/* General window events ----------------------------------------------- */
	public void windowClosed(WindowEvent e) {
	}

	public void windowDeactivated(WindowEvent e) {
	}

	public void focusLost(FocusEvent e) {
	}

	public void windowDeiconified(WindowEvent e) {
	}

	public void windowIconified(WindowEvent e) {
	}

	public void windowOpened(WindowEvent e) {
	}

	public void windowActivated(WindowEvent e) {
	}

	/**
	 * housekeepin'
	 * 
	 * @param e
	 *            windowclosing event (unused)
	 */
	public void windowClosing(WindowEvent e) {
		if (closed)
			return;
		if (isChangeSaved()) {
			Xmipp_Tomo.ExitValues choice = dialogYesNoCancel("Close window",
					"Are you sure you want to close this window?");
			switch (choice) {
			case NO:
			case CANCEL:
				return;
			}
		} else {
			Xmipp_Tomo.ExitValues choice = dialogYesNoCancel("Save changes",
					"Do you want to save changes?");
			switch (choice) {
			case YES:
				String path = dialogSave();
				if ("".equals(path))
					return;
				saveFile(getModel(), path);
				break;
			case CANCEL:
				return;
			}
		}

		// close the window (no return back...)
		setVisible(false);
		setModel(null);
		dispose();
		WindowManager.removeWindow(this);
		closed = true;
	}

	/* Mouse events ----------------------------------------------- */

	// based on ij.gui.ImageCanvas
	public void mouseMoved(MouseEvent e) {
		// update status when moving more that 12 pixels
		// if(getCursorDistance(e.getX(),e.getY())> 144)
		updateStatusText();
		setCursorLocation(e.getX(), e.getY());

	}

	public void mouseExited(MouseEvent e) {
	}

	public void mouseDragged(MouseEvent e) {
	}

	public void mouseClicked(MouseEvent e) {
	}

	public void mouseEntered(MouseEvent e) {
	}

	public void mouseReleased(MouseEvent e) {
		Xmipp_Tomo.debug("Mouse released");
	}

	public void mousePressed(MouseEvent e) {
	}

	// Property changes
	// Handle changes in the number of projections (either first load, or reload
	// prior to saving)
	public void propertyChange(PropertyChangeEvent event) {
		// Xmipp_Tomo.debug(event.getNewValue().toString());
		if (TomoData.Properties.NUMBER_OF_PROJECTIONS.name().equals(
				event.getPropertyName()))
			if (isReloadingFile())
				setStatus("Recovering original file... "
						+ getStatusChar((Integer) (event.getNewValue())));
			else
				setNumberOfProjections((Integer) (event.getNewValue()));
	}

	/* Getters/ setters --------------------------------------------------- */

	private ImageIcon getIcon(String name) {
		ImageIcon res = icons.get(name);
		if (res == null) {
			java.net.URL imageURL = Xmipp_Tomo.class.getResource(name);
			if (imageURL != null) {
				res = new ImageIcon(imageURL);
				icons.put(name, res);
			}
		}
		return res;
	}

	/**
	 * @return the model
	 */
	public TomoData getModel() {
		return model;
	}

	/**
	 * @param model
	 *            the model to set
	 */
	public void setModel(TomoData model) {
		this.model = model;

		if (this.model != null)
			// model change events
			this.model.addPropertyChangeListener(this);

	}

	private JLabel getStatusLabel() {
		return statusLabel;
	}

	private char getStatusChar(int i) {
		String animation = "\\|/-\\|/-";
		return animation.charAt(i % animation.length());
	}

	void setCursorLocation(int x, int y) {
		cursorLocation.setLocation(x, y);
	}

	public Rectangle getMaximumBounds() {
		return maxWindowSize;

	}

	// no need to reduce insets size (like ImageWindow does) - try to use
	// Container.insets (instead of ImageWindow.getInsets())
	public Insets getInsets() {
		// Insets insets = new Insets(0,0,0,0);

		return insets();
	}

	/**
	 * Update all window components related to the number of projections. The
	 * number itself is stored in the model
	 * 
	 * @param n
	 *            - see Tomodata.numberOfProjections
	 */
	public void setNumberOfProjections(int n) {
		if (projectionScrollbar != null)
			projectionScrollbar.setMaximum(n);
	}

	private int getCursorX() {
		return cursorLocation.x;
	}

	private int getCursorY() {
		return cursorLocation.y;
	}

	private double getCursorValue() {
		return getModel().getPixelValue(getCursorX(), getCursorY());
	}

	/*
	 * It requires that the model has been set (with setModel() )
	 */
	public String getCursorValueAsString() {
		if (getModel() == null)
			return "";

		String cursorValue = String.valueOf(getCursorValue());
		// include the decimal point in the count
		int limit = Math.min(cursorValue.length(), MAX_DIGITS + 1);
		return cursorValue.substring(0, limit);
	}

	public String getTitle() {
		String title = TITLE;
		if (getModel() != null)
			title = title + " > " + getModel().getFileName() + " ("
					+ getModel().getWidth() + "x" + getModel().getHeight()
					+ ")";
		return title;
	}

	/**
	 * Show a file open browser and...
	 * 
	 * @return the path of the file chosen by the user, or empty string (not
	 *         null) if the user cancels the dialog
	 * @deprecated Use {@link FileDialog#dialogOpen()} instead
	 */
	private String dialogOpen() {
		return FileDialog.openDialog("Import...", this);
	}

	/**
	 * Show a file save browser and...
	 * 
	 * @return the path of the file chosen by the user, or empty string (not
	 *         null) if the user cancels the dialog
	 */
	public String dialogSave() {
		return XrayImportDialog.dialogSave(getModel().getFilePath(), ".mrc");
	}

	/**
	 * @param title
	 *            Dialog title
	 * @param message
	 *            Dialog message
	 * @return Xmipp_Tomo.ExitValues.CANCEL / YES / NO
	 */
	private Xmipp_Tomo.ExitValues dialogYesNoCancel(String title, String message) {
		return XrayImportDialog.dialogYesNoCancel(title, message);
	}

	public static void alert(String message) {
		IJ.error("XmippTomo", message);
	}

	public void saveFile(TomoData model, String path) {
		setStatus("Saving...");
		model.setFile(path);
		new TiltSeriesIO(false).write(model);
		setStatus("Done");
		setChangeSaved(true);
	}

	/********************** Main - class testing **************************/

	private void gridBagLayoutTest() {
		/*
		 * Container pane=getContentPane();
		 * 
		 * JButton button; pane.setLayout(new GridBagLayout());
		 * GridBagConstraints c = new GridBagConstraints();
		 * 
		 * button = new JButton("Long-Named Button 4"); c.fill =
		 * GridBagConstraints.HORIZONTAL; c.ipady = 40; //make this component
		 * tall c.weightx = 0.0; c.gridwidth = 3; c.gridx = 0; c.gridy = 0;
		 * pane.add(button, c);
		 * 
		 * setVisible(true);
		 */
	}

	private void layoutTest() {
		addMainPanels();
		setVisible(true);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		TomoWindow w = new TomoWindow();
		/*
		 * frame.setSize(400, 400); frame.setLocationRelativeTo(null);
		 * frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		 */

		// w.gridBagLayoutTest();
		w.layoutTest();

	}

	/**
	 * @return the firstAction
	 */
	public UserAction getFirstAction() {
		return firstAction;
	}

	/**
	 * @param firstAction
	 *            the firstAction to set
	 */
	public void setFirstAction(UserAction firstAction) {
		this.firstAction = firstAction;
	}

	/**
	 * @return the lastAction
	 */
	public UserAction getLastAction() {
		if (lastAction == null)
			lastAction = firstAction;
		return lastAction;
	}

	/**
	 * @param lastAction
	 *            the lastAction to set
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
	 * @param windowId
	 *            the windowId to set
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
	 * @param plugin
	 *            the plugin to set
	 */
	public void setPlugin(Plugin plugin) {
		this.plugin = plugin;
	}

	/**
	 * @deprecated
	 * @return the reloadingFile
	 */
	private boolean isReloadingFile() {
		return (getLastCommandState() == Command.State.RELOADING);
	}

	/**
	 * @deprecated
	 * @param reloadingFile
	 *            the reloadingFile to set
	 */
	public void setReloadingFile(boolean reloadingFile) {
		if(reloadingFile)
			setLastCommandState(Command.State.RELOADING);
		else
			setLastCommandState(Command.State.LOADED);
	}

	/**
	 * @return the changeSaved
	 */
	private boolean isChangeSaved() {
		return changeSaved;
	}

	/**
	 * @param changeSaved
	 *            the changeSaved to set
	 */
	public void setChangeSaved(boolean changeSaved) {
		this.changeSaved = changeSaved;
	}

	/**
	 * @return the playing
	 */
	public boolean isPlaying() {
		return playing;
	}

	/**
	 * @param playing
	 *            the playing to set
	 */
	public void setPlaying(boolean playing) {
		this.playing = playing;
	}

	/*
	 * public void paint(Graphics g) { paintComponents(g); // paintAll(g);
	 * /*menuPanel.paint(g); viewPanel.paint(g); statusPanel.paint(g); //
	 * root.paint(g); }
	 */

	/**
	 * @return the canvas
	 * 
	 *         public ImageCanvas getCanvas() { return canvas; }
	 */

	/**
	 * @param canvas
	 *            the canvas to set
	 */
	private void setCanvas(ImageCanvas canvas) {
		// this.canvas = canvas;
		ic = canvas;
	}

	/**
	 * @return the timer
	 */
	private Timer getTimer() {
		if (timer == null) {
			timer = new Timer(2000, this);
			timer.setRepeats(false);
		}
		return timer;
	}

	/**
	 * @param timer
	 *            the timer to set
	 */
	private void setTimer(Timer timer) {
		this.timer = timer;
	}

	/**
	 * @return the projectionScrollbar
	 */
	public LabelScrollbar getProjectionScrollbar() {
		return projectionScrollbar;
	}

	/**
	 * @param projectionScrollbar
	 *            the projectionScrollbar to set
	 */
	public void setProjectionScrollbar(LabelScrollbar projectionScrollbar) {
		this.projectionScrollbar = projectionScrollbar;
	}

	/**
	 * @return the controller
	 */
	public TomoController getController() {
		return controller;
	}

	/**
	 * @param controller
	 *            the controller to set
	 */
	public void setController(TomoController controller) {
		this.controller = controller;
	}

	public void setLastCommandState(Command.State lastCommandState) {
		this.lastCommandState = lastCommandState;
	}

	public Command.State getLastCommandState() {
		return lastCommandState;
	}
}