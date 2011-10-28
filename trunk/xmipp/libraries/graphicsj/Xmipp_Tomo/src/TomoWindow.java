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
import ij.gui.Roi;
import java.awt.*;
import java.awt.event.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.Hashtable;
import java.util.LinkedList;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.Timer;
import javax.swing.UIManager;
import javax.swing.plaf.FontUIResource;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeSelectionModel;

/**
 * - Why?
 * MVC paradigm -> you need a View. Hence TomoWindow
 * 
 * TODO: -low- refactor Swing code to a superclass
 * 
 * implements: 
 * - WindowListener: windowclosing and the like
 * - MouseMotionListener: cursor position
 * - PropertyChangeListener: model properties changes
 * - DialogListener to capture GenericDialogs,
 * - ComponentListener to handle window resizing
 * - ActionListener: for timers (@see getTimer())
 */
public class TomoWindow extends ImageWindow implements WindowListener,
		MouseMotionListener, ActionListener,
		PropertyChangeListener, DialogListener, ComponentListener {
	/**
	 * Protocol for adding new buttons to the workflow/UI:
	 * - declare its command in Command
	 * - define the action method in TomoController
	 * - include the command in its menu list (for example, commandsMenuFile) 
	 * - set a proper enable/disable cycle along user interaction
	 * (tipically start disabled, enable it when user can click it)
	 */
	
	// TODO: monitor progress of commands stdout / log file inside workflowview / etc.
	// 		Implement with SwingWorker? (@see example in Swingworker documentation)
	// TODO: font antialiasing - @see ZZZZ
	// TODO: zoom, scroll... (quite buggy right now...)
	// TODO: remove load canceling (in the new approach with a cache, slices will be loaded one by one)
	// TODO: preprocessing - non-blocking previews, i.e., allow the user to change the current slice with the preproc dialog open
	// TODO: preproc.export or option in save dialog - write a file without the disabled (discarded) projections, not simply setting enable in selfile

	// for serialization only - just in case
	private static final long serialVersionUID = -4063711975454855701L;

	// MVC Controller. Controller Should not be needed here (passive window)
	private TomoController controller;

	
	// total number of digits to display as a String
	private static int MAX_DIGITS = 9;
	// Window title (without extra data; see getTitle() )
	private final static String TITLE = "XmippTomo";
	// minimum dimension of tabbed menu panel (window top), needed to adjust automatically the window size to
	// the tabbed menu needs (so it does not turn to a scrollable tabbed menu)
	private static int MENUPANEL_MINWIDTH = 600, MENUPANEL_MINHEIGHT = 150;
	// maximum window size
	private static Rectangle maxWindowSize = new Rectangle(800, 800);
	// miliseconds to wait for a Plugin dialog to display (so the dialog capture
	// can start)
	public static int WAIT_FOR_DIALOG_SHOWS = 800;

	public static String PLAY_ICON = "resources/icon-play.png",
			PAUSE_ICON = "resources/icon-pause.png";
	
	private Hashtable<String, ImageIcon> icons = new Hashtable<String, ImageIcon>();

	private Hashtable<String, JButton> buttons = new Hashtable<String, JButton>();
	private Hashtable<String, JCheckBox> checkboxes = new Hashtable<String, JCheckBox>();

	// Lists of commands for each menu tab
	private static java.util.List<Command> commandsMenuFile = new LinkedList<Command>() {
		{
			add(XmippTomoCommands.LOAD);
			add(XmippTomoCommands.XRAY);
			add(XmippTomoCommands.CONVERT);
			add(XmippTomoCommands.DEFINE_TILT);
			add(XmippTomoCommands.DISCARD_PROJECTION);
		}
	};

	private static java.util.List<Command> commandsMenuPreproc = new LinkedList<Command>() {
		{
			add(XmippTomoCommands.NORMALIZE_SERIES);
			add(XmippTomoCommands.HOTSPOT_REMOVAL);
			add(XmippTomoCommands.GAUSSIAN);
			add(XmippTomoCommands.MEDIAN);
			add(XmippTomoCommands.SUB_BACKGROUND);
			add(XmippTomoCommands.ENHANCE_CONTRAST);
			add(XmippTomoCommands.GAMMA_CORRECTION);
			add(XmippTomoCommands.BANDPASS);
			add(XmippTomoCommands.HISTOGRAM_EQUALIZATION);
			add(XmippTomoCommands.CROP);
			// add(XmippTomoCommands.APPLY);
		}
	};

	private static java.util.List<Command> commandsMenuAlign = new LinkedList<Command>() {
		{
			add(XmippTomoCommands.ALIGN_AUTO);
			add(XmippTomoCommands.ALIGN_MANUAL);
			add(XmippTomoCommands.ALIGN_CORRELATION);
		}
	};

	private static java.util.List<Command> commandsMenuDebug = new LinkedList<Command>() {
		{
			add(XmippTomoCommands.PRINT_WORKFLOW);
			add(XmippTomoCommands.CURRENT_PROJECTION_INFO);
		}
	};

	// which commands to enable after loading a file
	private static java.util.List<Command> commandsEnabledAfterLoad = new LinkedList<Command>() {
		{
			// disabled until native writing (using Xmipp library) is implemented
			//add(Command.SAVE);
			add(XmippTomoCommands.NORMALIZE_SERIES);
			add(XmippTomoCommands.DEFINE_TILT);
			add(XmippTomoCommands.DISCARD_PROJECTION);
			add(XmippTomoCommands.GAUSSIAN);
			add(XmippTomoCommands.MEDIAN);
			add(XmippTomoCommands.SUB_BACKGROUND);
			add(XmippTomoCommands.ENHANCE_CONTRAST);
			add(XmippTomoCommands.GAMMA_CORRECTION);
			add(XmippTomoCommands.BANDPASS);
			add(XmippTomoCommands.HISTOGRAM_EQUALIZATION);
			add(XmippTomoCommands.CROP);
			//add(Command.HOTSPOT_REMOVAL);
			// TODO: convert needs to be adapted to the new image-by-image I/O model
			// add(XmippTomoCommands.CONVERT);
			add(XmippTomoCommands.PRINT_WORKFLOW);
			add(XmippTomoCommands.CURRENT_PROJECTION_INFO);
			add(XmippTomoCommands.MEASURE);
			add(XmippTomoCommands.PLAY);
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
	// TODO: not needed - remove changeSaved
	// true if changes are saved (so you can close the window without pain)
	private boolean changeSaved = true;
	
	// window identifier - useful when you've got more than 1 window
	private int windowId = -1;

	/* Window GUI components */
	/* 4 main panels */
	private JTabbedPane menuPanel;
	private JPanel imagePanel, stackPanel, viewControlsPanel,controlPanel, statusPanel;
	private JPanel viewsPanel;

	private WorkflowView projectView; // alternative tree-like views: Prefuse
	
	// hack for reusing ImageJ ImageWindow
	// TODO: verify if the hack is still needed
	private JFrame realWindow;

	// Text fields and labels
	private JTextField tiltTextField;
	private JLabel statusLabel;

	// Control scrollbars
	private LabelScrollbar projectionScrollbar;



	Point cursorLocation = new Point();

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

		realWindow.addWindowListener(this);
		realWindow.addComponentListener(this);
		// EXIT_ON_CLOSE finishes ImageJ too...
		// DO NOTHING allows for closing confirmation dialogs and the like
		realWindow.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		WindowManager.addWindow(this);
		realWindow.setTitle(getTitle());
		addMainPanels(getController().getWorkflow());
	}

	public TomoWindow(int windowId) {
		this();
		setWindowId(windowId);
	}

	// Interface to the JFrame we are really displaying to
	private Container getContentPane() {
		return realWindow.getContentPane();
	}

	private void addButton(Command cmd, Container panel) throws Exception {
		Action a = new Action(getController(), cmd);
		JButton b = new JButton(a);
		panel.add(b);
		buttons.put(cmd.getId(), b);
	}
	
	private JCheckBox addCheckbox(Command cmd, Container panel) throws Exception {
		Action a = new Action(getController(), cmd);
		JCheckBox b = new JCheckBox(a);
		b.setLabel(cmd.getLabel());
		panel.add(b);
		checkboxes.put(cmd.getId(), b);
		return b;
	}
	
	private void addCheckbox(Command cmd, Container panel,boolean checked) throws Exception {
		addCheckbox(cmd, panel).setSelected(checked);
	}

	public JButton getButton(String id) {
		return buttons.get(id);
	}
	
	public JCheckBox getCheckbox(String id) {
		return checkboxes.get(id);
	}


	private void addMainPanels(Workflow workflow) {
		realWindow.setLayout(new BoxLayout(getContentPane(),
				BoxLayout.PAGE_AXIS));

		// BoxLayout component order is sequential, so all the panels
		// MUST be added in order (even if they are empty)

		// TABBED MENU PANEL
		menuPanel = new JTabbedPane();
		try{
			addMenuTabs();
		}catch (Exception ex){
			Logger.debug("Tomowindow.addMainPanels - addmenutabs ", ex);
		}
		menuPanel.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
		menuPanel.setPreferredSize(new Dimension(MENUPANEL_MINWIDTH,
				MENUPANEL_MINHEIGHT));
		getContentPane().add(menuPanel);

		// VIEW & CONTROLS PANELS
		viewsPanel = new JPanel();
		viewsPanel.setLayout(new BorderLayout());
		viewControlsPanel = new JPanel();
		stackPanel = new JPanel();
		stackPanel.setLayout(new BoxLayout(stackPanel,BoxLayout.PAGE_AXIS));
		imagePanel = new JPanel();
        stackPanel.setBorder(BorderFactory.createCompoundBorder(
        		BorderFactory.createTitledBorder("Stack View"),
                BorderFactory.createEmptyBorder(5,5,5,5)));
        stackPanel.add(imagePanel);
        stackPanel.add(viewControlsPanel);
		viewsPanel.add(stackPanel,BorderLayout.CENTER);
		
		projectView=new WorkflowView(workflow,null);
		projectView.addPropertyChangeListener(this);
		
		viewsPanel.add(projectView,BorderLayout.EAST);
		getContentPane().add(viewsPanel);
		
		controlPanel = new JPanel();
		getContentPane().add(controlPanel);

	
		if (getStackModel() != null) {
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
		if (getStackModel() == null)
			return;

		// remove previous canvas if any
		for (int i = 0; i < imagePanel.getComponentCount(); i++)
			imagePanel.remove(i);

		setCanvas(new TomoImageCanvas(getStackModel()));
		
		imagePanel.setLayout(new ImageLayout(getCanvas()));
		imagePanel.add(getCanvas());
		getCanvas().addMouseMotionListener(this);
		// listen to ImageWindow properties changes
		// TODO: addView - is really needed?
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
			projectionScrollbar.setMaximum(getStackModel().getNumberOfProjections());
			return;
		}
		
		/* remove previous controls if any
		for (int i = 0; i < controlPanel.getComponentCount(); i++)
			controlPanel.remove(i); */

		FlowLayout fl2 = new FlowLayout();
		fl2.setAlignment(FlowLayout.CENTER);
		viewControlsPanel.setLayout(fl2);
		controlPanel.setLayout(fl2);

		JLabel tiltTextLabel = new JLabel("Tilt");
		tiltTextField = new JTextField(5);
		// the user won't change angles one by one - he will change all of them
		// at once with the tilt change dialog
		tiltTextField.setEditable(false);
		viewControlsPanel.add(tiltTextLabel);
		viewControlsPanel.add(tiltTextField);

		projectionScrollbar = new LabelScrollbar(1, getStackModel()
				.getNumberOfProjections());
		projectionScrollbar.setText("Projection #");
		projectionScrollbar.addAdjustmentListener(getController());
		viewControlsPanel.add(projectionScrollbar);

		try {
			addButton(XmippTomoCommands.PLAY, viewControlsPanel);
			addCheckbox(XmippTomoCommands.PLAY_LOOP, viewControlsPanel,true);
			addButton(XmippTomoCommands.MEASURE, controlPanel);
			addButton(XmippTomoCommands.ADJUSTBC, controlPanel);
		} catch (Exception ex) {
			Logger.debug("addControls - play");
		}
		realWindow.pack(); // adjust window size to host the new controls panel

	}

	/**
	 * Add menu tabs populated with their buttons (specified in the lists at the
	 * beginning of the class)
	 */
	private void addMenuTabs() throws Exception {
		// File
		JPanel fileTabPanel = new JPanel(false); // false = no double buffer (flicker)
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

	// TODO: (bug) sometimes the label does not refresh automatically (for example, when saving a file)
	// probably was caused by not running some processes in background... check if this bug still happens
	public void setStatus(String text) {
		getStatusLabel().setText(text);
	}

	private void updateStatusText() {
		// current/total projections are shown in the projection scrollbar
		// itself
		if (getStackModel() != null)
			setStatus("x = " + getCursorX() + ", y = " + getCursorY()
					+ ", value = "+ getCursorValueAsString()); 
	}

	int getCursorDistance(int x, int y) {
		return (int) cursorLocation.distance(x, y);
	}

	// TODO: -current- zoom out big images
	public void refreshImageCanvas() {
		if (getCanvas() != null) {	
			ImagePlus img=getStackModel().getCurrentImage();
			if(img != null){
				setImage(img);
				// required for IJ actions to work
				WindowManager.setTempCurrentImage(img);
				// scale image if needed
				 if(img.getWidth() > 512 || img.getHeight() > 512){
					// TODO: take into consideration height too (for non-square images)
					double mag=(512.0/img.getWidth()) * 100;
					ZoomPlugin plugin = new ZoomPlugin(mag);
					plugin.run(img);
				}
			}
			// TODO: make the image canvas size fixed - next line is not working for that
			getCanvas().setPreferredSize(new Dimension(512,512));
			getCanvas().setImageUpdated();
			getCanvas().repaint();
			

		}
	}

	public void enableButtonsAfterLoad() {
		for (Command c : commandsEnabledAfterLoad) {
			getButton(c.getId()).setEnabled(true);
		}
	}

	public void setImagePlusWindow() {
		imp = getStackModel().getCurrentImage();
		imp.setWindow(this);
	}

	/*
	 * Timer events
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
/**
 * @deprecated
 */
 	public void addUserAction(UserAction action) {
 
		projectView.newOperation(action);
	}

	public void changeIcon(String buttonId, String iconName) {
		getButton(buttonId).setIcon(getIcon(iconName));
	}
	
	public void changeLabel(String buttonId, String label) {
		getButton(buttonId).setText(label);
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
//		Component c = e.getComponent();
//		// Xmipp_Tomo.debug("componentResized event from " +
//		// c.getClass().getName());
//		// + "; new size: " + c.getSize().width + ", "
//		// + c.getSize().height + " - view size:" + viewPanel.getSize());
//		getTimer().stop();
//		// getTimer().setDelay(2000);
//		timer.start();
//		resizeView();
	}

	/**
	 * Fit the image canvas to the space available in the view, keeping the
	 * aspect ratio
	 */
//	public void resizeView() {
//		// TODO: resizeView does not work... - for window resizing, call pack() with Timers?
//		if ((imagePanel == null) || (getCanvas() == null))
//			return;
//
//		double factorWidth = imagePanel.getWidth()
//				/ getCanvas().getPreferredSize().getWidth();
//		double factorHeight = imagePanel.getHeight()
//				/ getCanvas().getPreferredSize().getHeight();
//		double factor = 1;
//
//		if (factorWidth < 1) {
//			if (factorHeight < 1)
//				// apply maximum of both (since both are < 1 )
//				factor = Math.max(factorWidth, factorHeight);
//			else
//				factor = factorWidth;
//		} else {
//			if (factorHeight < 1)
//				factor = factorHeight;
//			else
//				factor = Math.min(factorWidth, factorHeight);
//		}
//
//		// if size change is minimal or none, don't resize canvas
//		if (Math.abs(factor - 1) < 0.03)
//			return;
//
//		// Xmipp_Tomo.debug("resize: " + factorWidth + "," + factorHeight + ", "
//		// + factor );
//		resizeView(factor);
//
//	}
//
//	public void resizeView(double factor) {
//		double w = getCanvas().getWidth() * factor;
//		double h = getCanvas().getHeight() * factor;
//
//		// getCanvas().resizeCanvas((int)w, (int)h);
//		/*
//		 * ImageCanvas canvas = getModel().getImage().getCanvas();
//		 * getCanvas().setMagnification(factor); canvas.setSourceRect(new
//		 * Rectangle(0, 0, (int)(w/factor), (int)(h/factor)));
//		 * canvas.setDrawingSize((int)w, (int)h); realWindow.pack();
//		 * canvas.repaint();
//		 */
//		WindowManager.setTempCurrentImage(getStackModel().getCurrentImage());
//		// Xmipp_Tomo.debug("Set... " + "zoom="+ ((int) (factor * 100)));
//		IJ.run("Set... ", "zoom=" + ((int) (factor * 100)));
//		refreshImageCanvas();
//
//	}

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
		// TODO: windowClosing - stop running threads (file load, for example) ??
		// TODO: now that each step is automatically written to disk, it makes not a lot of sense to ask the user
		// if he wants to save changes...
		if (closed)
			return;
		if (isChangeSaved()) {
			ExitValue choice = dialogYesNoCancel("Close window",
					"Are you sure you want to close this window?",false);
			switch (choice) {
				case NO:
				case CANCEL:
					return;
			}
		} else {
			ExitValue choice = dialogYesNoCancel("Save changes",
					"Do you want to save changes?",true);
			switch (choice) {
			case YES:
				String path = TomoFileDialog.saveDialog("Save...", this);
				if ("".equals(path))
					return;
				controller.saveFile(this, getModel(), path);
				break;
			case CANCEL:
				return;
			}
		}

		// close the window (no return back...)
		realWindow.setVisible(false);
		realWindow.dispose();
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

	// TODO: remove mouse handling events if not needed
	public void mouseExited(MouseEvent e) {
	}

	public void mouseDragged(MouseEvent e) {
	}

	public void mouseClicked(MouseEvent e) {
	}

	public void mouseEntered(MouseEvent e) {
	}

	public void mouseReleased(MouseEvent e) {
		Logger.debug("Mouse released");
	}

	public void mousePressed(MouseEvent e) {
	}

	// MVC - handle changes in the Model properties
	public void propertyChange(PropertyChangeEvent event) {
		// Xmipp_Tomo.debug(event.getNewValue().toString());
		if (TomoData.Properties.NUMBER_OF_PROJECTIONS.name().equals(event.getPropertyName())){
			if (isReloadingFile())
				setStatus("Recovering original file... "
						+ getStatusChar((Integer) (event.getNewValue())));
			else
				setNumberOfProjections((Integer) (event.getNewValue()));
		}else if (WorkflowView.ACTION_SELECTED.equals(event.getPropertyName()) ||
				StackModel.Properties.CURRENT_PROJECTION.name().equals(event.getPropertyName())){
			refreshImageCanvas();
			updateTitle();
		}else if(StackModel.Properties.CURRENT_PROJECTION_NUMBER.name().equals(event.getPropertyName())){
			updateCurrentTiltAngleText();
			getProjectionScrollbar().setValue(getStackModel().getCurrentProjectionNumber());
			refreshImageCanvas();	
			updateStatusText();
			// Discard button label
			if(getStackModel().isCurrentEnabled())
				changeLabel(XmippTomoCommands.DISCARD_PROJECTION.getId(), XmippTomoCommands.DISCARD_PROJECTION.getLabel());
			else
				changeLabel(XmippTomoCommands.DISCARD_PROJECTION.getId(), XmippTomoCommands.UNDO_DISCARD_PROJECTION.getLabel());
		}else if(TomoData.Properties.CURRENT_TILT_ANGLE.name().equals(event.getPropertyName())){
			updateCurrentTiltAngleText();
		}else if(TomoData.Properties.CURRENT_PROJECTION_ENABLED.name().equals(event.getPropertyName())){
		
			if(getStackModel().isCurrentEnabled())
				changeLabel(XmippTomoCommands.DISCARD_PROJECTION.getId(), XmippTomoCommands.DISCARD_PROJECTION.getLabel());
			else
				changeLabel(XmippTomoCommands.DISCARD_PROJECTION.getId(), XmippTomoCommands.UNDO_DISCARD_PROJECTION.getLabel());
			refreshImageCanvas();	
		}
	}

	private void updateCurrentTiltAngleText(){
		String text="";
		Double tilt=getStackModel().getCurrentTiltAngle();
		if(tilt != null)
			text=String.valueOf(tilt);
		setCurrentTiltAngleText(text);
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

	private TomoData getModel() {
		return getController().getModel();
	}

	private StackModel getStackModel() {
		return getController().getStackModel();
	}

	private JLabel getStatusLabel() {
		return statusLabel;
	}

	private char getStatusChar(int i) {
		String animation = "\\|/-\\|/-";
		return animation.charAt(i % animation.length());
	}

	
	private void setCurrentTiltAngleText(String tilt){
		if((tilt != null) && (tiltTextField != null))
			tiltTextField.setText(tilt);
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
	 * Update all window components related to the number of projections.  
	 * @param n
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
		return getStackModel().getPixelValue(getCursorX(), getCursorY());
	}

	/*
	 * It requires that the stack model has been set 
	 */
	public String getCursorValueAsString() {
		if (getStackModel() == null)
			return "";

		String cursorValue = String.valueOf(getCursorValue());
		// include the decimal point in the count
		int limit = Math.min(cursorValue.length(), MAX_DIGITS + 1);
		return cursorValue.substring(0, limit);
	}

	
	public String getTitle() {
		String title = TITLE;
		if (getStackModel() != null){
			String filename = getStackModel().getCurrentFileName();
			if(filename != null)
				title = title + " > " + filename;
			if(getStackModel().getCurrentWidth() != -1)
				title = title + " ("+ getStackModel().getCurrentWidth() + "x" + getStackModel().getCurrentHeight()
						+ ")," + getStackModel().getCurrentBitDepth() + "bits";
		}
		return title;
	}
	
	public void setTitle(String title){
		realWindow.setTitle(title);
	}
	
	public void updateTitle(){
		setTitle(getTitle());
	}

	/**
	 * @param title
	 *            Dialog title
	 * @param message
	 *            Dialog message
	 * @return Xmipp_Tomo.ExitValues.CANCEL / YES / NO
	 */
	private ExitValue dialogYesNoCancel(String title, String message,boolean showCancelButton) {
		return XrayImportDialog.dialogYesNoCancel(title, message,showCancelButton);
	}

	public static void alert(String message) {
		IJ.error("XmippTomo", message);
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

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		/*
		 * frame.setSize(400, 400); frame.setLocationRelativeTo(null);
		 * frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		 */

		// w.gridBagLayoutTest();

	}

	public int getWindowId() {
		return windowId;
	}
	


	public void setWindowId(int windowId) {
		this.windowId = windowId;
	}

	// TODO: move to background method as a local variable
	public Plugin getPlugin() {
		return plugin;
	}

	public void setPlugin(Plugin plugin) {
		this.plugin = plugin;
	}

	/**
	 * @deprecated
	 */
	private boolean isReloadingFile() {
		return (getLastCommandState() == Command.State.RELOADING);
	}

	/**
	 * @deprecated
	 */
	public void setReloadingFile(boolean reloadingFile) {
		if(reloadingFile)
			setLastCommandState(Command.State.RELOADING);
		else
			setLastCommandState(Command.State.LOADED);
	}

	private boolean isChangeSaved() {
		return changeSaved;
	}

	public void setChangeSaved(boolean changeSaved) {
		this.changeSaved = changeSaved;
	}


	private Timer getTimer() {
		if (timer == null) {
			timer = new Timer(2000, this);
			timer.setRepeats(false);
		}
		return timer;
	}

	private void setTimer(Timer timer) {
		this.timer = timer;
	}

	public LabelScrollbar getProjectionScrollbar() {
		return projectionScrollbar;
	}


	public void setProjectionScrollbar(LabelScrollbar projectionScrollbar) {
		this.projectionScrollbar = projectionScrollbar;
	}


	public TomoController getController() {
		return controller;
	}


	public void setController(TomoController controller) {
		this.controller = controller;
	}

	public void setLastCommandState(Command.State lastCommandState) {
		this.lastCommandState = lastCommandState;
	}

	public Command.State getLastCommandState() {
		return lastCommandState;
	}


	
	/**
	 * Avoid the image popping out to a new Image Window
	 * Remember to call unprotect after performing operations
	 */
	public void protectWindow(){
		getImagePlus().setWindow(null);
	}
	
	public void unprotectWindow(){
		getImagePlus().setWindow(this);
	}
	
	/**
	 * @deprecated
	 * @return
	 */
	public boolean isLoadCanceled(){
		return(getLastCommandState() == Command.State.CANCELED);
	}
	
	/**
	 * ImageWindow handles the canvas directly (with no setter)
	 * @param canvas
	 */
	private void setCanvas(ImageCanvas canvas) {
		ic = canvas;
	}
	

}