/***************************************************************************
 * Authors:     Juanjo Vega
 * 				J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
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

package xmipp.viewer.windows;

import ij.ImagePlus;

import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.swing.AbstractAction;
import javax.swing.AbstractButton;
import javax.swing.ActionMap;
import javax.swing.ButtonGroup;
import javax.swing.ComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;
import javax.swing.KeyStroke;
import javax.swing.LookAndFeel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.event.ListDataListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.JTableHeader;
import xmipp.ij.commons.ImagePlusLoader;
import xmipp.ij.commons.Tool;
import xmipp.ij.commons.XmippApplication;
import xmipp.ij.commons.XmippIJUtil;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.jni.MDRow;
import xmipp.utils.DEBUG;
import xmipp.utils.Param;
import xmipp.utils.QuickHelpJDialog;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippLabel;
import xmipp.utils.XmippMenuBarCreator;
import xmipp.utils.XmippMenuCreator;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.utils.XmippQuestionDialog;
import xmipp.utils.XmippResource;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.RowHeaderRenderer;
import xmipp.viewer.ctf.TasksEngine;
import xmipp.viewer.ctf.iCTFGUI;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.models.GalleryData;
import xmipp.viewer.models.GalleryRowHeaderModel;
import xmipp.viewer.models.ImageGalleryTableModel;
import xmipp.viewer.models.MetadataGalleryTableModel;
import xmipp.viewer.models.MicrographsTableModel;
import xmipp.viewer.particlepicker.extract.ExtractParticlePicker;
import xmipp.viewer.particlepicker.extract.ExtractPickerJFrame;
import xmipp.viewer.windows.ClassesJDialog;


/**
 * This is the main frame used in showj
 * @author airen
 *
 */
public class GalleryJFrame extends JFrame implements iCTFGUI
{
	private static final long serialVersionUID = -8957336972082018823L;

	private final static int DELAY_TO_UPDATE = 500;
	private static int update_counter = 0;
	// The following counter will be used to keep track of how many
	// windows are opened, the last one, should do System.exit
	// private static short windows_counter = 0;
	public ImageGalleryTableModel gallery;
	private GalleryRowHeaderModel rowHeaderModel;
	private int previousSelectedRow, previousSelectedCol;
	private JList rowHeader;
	// this flag will be used to avoid firing properties change events
	// when the change is from our code and not external user interaction
	private boolean isUpdating;
	private boolean autoAdjustColumns = false;
	private GalleryPopupMenu jpopUpMenuTable;
	private GalleryMenu menu;
	private XmippFileChooser fc;
	private SaveJDialog dlgSave = null;
	private boolean saved = false;
	private ClassesJDialog dlgClasses = null;

	private JLabel jlZoom;
	private JLabel jlGoToImage;
	private JLabel jlRows;
	private JLabel jlColumns;
	private JToggleButton jcbAutoAdjustColumns;
	private JButton btnChangeView;
	private JCheckBox jcbShowLabels;
	protected JPanel jpBottom;
	protected JSpinner jsColumns;
	protected JSpinner jsGoToImage;
	private JScrollPane jspContent;
	protected JSpinner jsRows;
	protected JSpinner jsZoom;
	// Components for combos
	protected JPanel cbPanel;
	protected JComboBox jcbBlocks;
	protected JComboBox jcbVolumes;
	protected JLabel jlBlocks;
	protected JLabel jlVolumes;
	protected TasksEngine ctfTasks;
	// private javax.swing.JToggleButton jtbNormalize;
	// private javax.swing.JToggleButton jtbUseGeometry;
	private JTable table;
	private JToolBar toolBar;
	private int width = -1;
	private JButton reslicebt;
	private String[] reslices = new String[] { "Z Negative (Front)", "Y Negative (Top)", "X Negative (Left)", "Y Positive (Bottom)",
			"X Positive (Right)" };;

	protected static final float MAX_HEIGHT_RATE = 2.0f / 3.0f;
	// this rate is width/height
	protected static final float DIM_RATE = 4.0f / 3.0f;
	protected static final int MIN_WIDTH = 650;
	protected static int MIN_HEIGHT;
	protected static int MAX_HEIGHT;
	protected static int MAX_WIDTH;
	protected static Dimension screenSize;
	private Integer rows, columns;
        protected JPanel commandspn;

	/** Store data about visualization */
	GalleryData data;
	private ExtractPickerJFrame extractframe;
	private ButtonGroup reslicegroup;

	private JComboBox imagecolumnscb;

	private JMenuItem rendercolumnmi;

	private Hashtable<String, ColumnInfo> imagecolumns;
	/** Some static initialization for fancy default dimensions */
	static
	{
		screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		float aux = (float) screenSize.height * MAX_HEIGHT_RATE;
		MAX_HEIGHT = Math.round(aux);
		aux = (float) MIN_WIDTH / DIM_RATE;
		MIN_HEIGHT = Math.round(aux);
		aux = (float) MAX_HEIGHT * DIM_RATE;
		MAX_WIDTH = Math.round(aux);
	}
    

	/** Initialization function after GalleryData structure is created */
	private void init(GalleryData data)
	{
		try
		{
			this.data = data;
			createModel();
			createGUI();
			XmippApplication.addInstance();
		}
		catch (Exception e)
		{
			DEBUG.printException(e);
		}
	}

	/** Constructors */
	public GalleryJFrame(String filename, Param parameters)
	{
		super();
		init(new GalleryData(this, filename, parameters, null));
	}

	public GalleryJFrame(String filename, MetaData md, Param parameters)
	{
		super();
		init(new GalleryData(this, filename, parameters, md));
	}

	/**
	 * Open another metadata separataly *
	 */
	public void openMetadata(MetaData md)
	{
		new GalleryJFrame(null, md, new Param());
	}

	/**
	 * Function to create the gallery type depending on the filename
	 * 
	 * @throws Exception
	 */
	private void createModel() throws Exception
	{
		gallery = data.createModel();
	}

	public GalleryData getData()
	{
		return data;
	}

	/** Close the application, check if changes first */
	public void close()
	{
		if (proceedWithChanges())
		{
			setVisible(false);
			dispose();
			XmippApplication.removeInstance();

		}
	}// function close

	/** Check if there are changes to proceed */
	public boolean proceedWithChanges()
	{
		boolean proceed = true;
		if (data.hasMdChanges())
		{
			XmippQuestionDialog dlg = new XmippQuestionDialog(GalleryJFrame.this, "Do you want to save metadata changes?");
			if (dlg.showDialog())
				try
				{
					save();
				}
				catch (Exception e)
				{
					showException(e);
				}
			else
				proceed = !dlg.isCanceled();
		}
		return proceed;
	}

	/** Set the title of the main windows depending on the gallery */
	private void setGalleryTitle()
	{
		setTitle(gallery.getTitle());
	}

	/**
	 * Function to create general GUI base on a TableModel. It will use helper
	 * functions to create different components of the GUI
	 */
	private void createGUI()
	{
		// Create file chooser and set current dir
		setIconImage(XmippResource.getIcon("xmipp_logo.png").getImage());
		if (data.getFileName() != null)
			fc = new XmippFileChooser(new File(data.getFileName()));
		else
			fc = new XmippFileChooser();
		ctfTasks = new TasksEngine(GalleryJFrame.this);

		isUpdating = true; // avoid handling some changes events

		setGalleryTitle();
		setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
		setMinimumSize(new Dimension(MIN_WIDTH, MIN_HEIGHT));
		addWindowListener(new WindowAdapter()
		{
			@Override
			public void windowClosing(WindowEvent arg0)
			{
				close();
			}
		});

		// Get main pane and set layout
		Container pane = getContentPane();
		JPanel container = new JPanel(new GridBagLayout());
		pane.add(container);
		// container.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();

		// Create toolbar buttons
		createToolbar();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridy = 0;
		container.add(toolBar, c);
		setInitialValues();

		// Create combos for selection of blocks and/or volumes
		createCombos();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridy = 1;
		container.add(cbPanel, c);
		updateVisibleCombos();

		jspContent = new GalleryScroll();
		// Create table
		createTable();
		c.fill = GridBagConstraints.BOTH;
		c.gridx = 0;
		c.gridy = 2;
		c.weightx = 1.0;
		c.weighty = 1.0;
		container.add(jspContent, c);
                
                commandspn = new JPanel(new FlowLayout(FlowLayout.RIGHT));
               
                c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridy = 3;
		c.weightx = 0;
		c.weighty = 0;
		container.add(commandspn, c);
				
		// Create the menu for table
		menu = new GalleryMenu();
		setJMenuBar(menu.getMenuBar());
		jpopUpMenuTable = new GalleryPopupMenu();
		menu.update();

		// pack();
		isUpdating = false;

		// Zoom in with Ctrl + P
		InputMap imap = container.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
		ActionMap amap = container.getActionMap();
		imap.put(KeyStroke.getKeyStroke("ctrl released P"), "zoomIn");
		amap.put("zoomIn", new AbstractAction()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				zoomChange(true);
			}

		});
		// Zoom in with Ctrl + O
		imap.put(KeyStroke.getKeyStroke("ctrl released M"), "zoomOut");
		amap.put("zoomOut", new AbstractAction()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				zoomChange(false);
			}

		});

		// Change view with Ctrl + Tab
		imap.put(KeyStroke.getKeyStroke("ctrl released I"), "changeView");
		amap.put("changeView", new AbstractAction()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				data.changeMode();
				reloadTableData();
			}
		});

		pack();
		XmippWindowUtil.centerWindows(this);
		setVisible(true);
		SwingUtilities.invokeLater(new Runnable()
		{
			public void run()
			{
				addComponentListener(new java.awt.event.ComponentAdapter()
				{
					public void componentResized(java.awt.event.ComponentEvent evt)
					{
						formComponentResized(evt);
					}
				});
			}
		});
	}

	private void setInitialValues()
	{
		boolean adjust = false;
		if (data.parameters.columns > 0)
			gallery.setColumns(data.parameters.columns);
		else if (data.parameters.rows > 0)
			gallery.setRows(data.parameters.rows);
		else if (!data.isRotSpectraMode())
			adjust = true;

		if (data.isMicrographsMode())
		{
			// setExtendedState(JFrame.MAXIMIZED_BOTH);
			width = screenSize.width - 50;
			int h = screenSize.height - 100;
			setPreferredSize(new Dimension(width, h));
		}
		else
		{
			int desiredCols = adjust ? (int) Math.ceil(Math.sqrt(gallery.getSize())) : gallery.getColumnCount();
			width = desiredCols * gallery.cellDim.width + 50;
			width = Math.min(Math.max(width, MIN_WIDTH), MAX_WIDTH);
			if (adjust)
			{
				gallery.adjustColumn(width - 50);
			}
			int h = gallery.getRowCount() * gallery.cellDim.height;
			h = Math.min(Math.max(h, MIN_HEIGHT), MAX_HEIGHT);
			setPreferredSize(new Dimension(width, h));
		}
		setAutoAdjustColumns(adjust);
	}

	/** Some tweaks over traditional JTable */
	public class GalleryScroll extends JScrollPane
	{
	}// class GalleryTable

	private void createTable()
	{
		// Create row header for enumerate rows
		try
		{
			rowHeaderModel = (data.md.isColumnFormat() || !data.isTableMode()) ? new GalleryRowHeaderModel(gallery.getRowCount(), 1)
					: new GalleryRowHeaderModel(data);
		}
		catch (Exception e1)
		{
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		rowHeader = new JList();
		rowHeader.setModel(rowHeaderModel);
		LookAndFeel.installColorsAndFont(rowHeader, "TableHeader.background", "TableHeader.foreground", "TableHeader.font");
		rowHeader.setCellRenderer(new RowHeaderRenderer());
		jspContent.setRowHeaderView(rowHeader);

		table = new JTable()
		{
			protected JTableHeader createDefaultTableHeader()
			{
				return gallery.getTableHeaderModel();
			}
		};
		// Create column model
		table.setColumnModel(gallery.getColumnModel());
		table.setModel(gallery);
		// int h = 25;
		// table.setRowHeight(h);
		// rowHeader.setFixedCellHeight(h);
		jspContent.setViewportView(table);

		gallery.setupTable(table);
		table.setRowSelectionAllowed(true);
		// DEBUG.printMessage("WIDTH: " + jspContent.getVisibleRect().width);
		// DEBUG.printMessage("preferred: " + getPreferredSize().toString());
		gallery.addTableModelListener(new TableModelListener()
		{

			@Override
			public void tableChanged(TableModelEvent e)
			{
				updateTable();
			}
		});

		table.addMouseListener(new java.awt.event.MouseAdapter()
		{
			public void mouseClicked(java.awt.event.MouseEvent evt)
			{
				tableMouseClicked(evt);
			}
		});

		// Zoom with Shift + MouseWeel
		table.addMouseWheelListener(new MouseWheelListener()
		{
			@Override
			public void mouseWheelMoved(MouseWheelEvent evt)
			{
				if (evt.isShiftDown())
					zoomChange(evt.getWheelRotation() < 0);
				else
					table.getParent().dispatchEvent(evt);
			}
		});

		// Add listener to recognize UP and DOWN keys
		table.addKeyListener(new KeyListener()
		{
			@Override
			public void keyTyped(KeyEvent arg0)
			{
			}

			@Override
			public void keyReleased(KeyEvent arg0)
			{
			}

			@Override
			public void keyPressed(KeyEvent arg0)
			{
				int vdir = 0, hdir = 0;

				switch (arg0.getKeyCode())
				{
				case KeyEvent.VK_DELETE:
					removeObjects(true);
					break;
				case KeyEvent.VK_UP:
					vdir = -1;
					break;
				case KeyEvent.VK_DOWN:
					vdir = 1;
					break;
				case KeyEvent.VK_LEFT:
					hdir = -1;
					break;
				case KeyEvent.VK_RIGHT:
					hdir = 1;
					break;
				}
				if (vdir != 0 || hdir != 0)
				{
					int newRow = table.getSelectedRow() + vdir;
					int col = table.getSelectedColumn() + hdir;
					//System.err.format("newRow: %d, col: %d\n", newRow, col);
					selectItem(newRow, col);

				}
			}// function keyPressed
		});

		updateViewState();

		if (!adjustColumns())
			updateTable(); // update table if columns have not changed
	}// function createTable

	private void zoomChange(boolean increase)
	{
		
		int deltha = increase ? 10 : -10;
		jsZoom.setValue((Integer) jsZoom.getValue() + deltha);
	}

	private void updateTable()
	{
		// if (table.isShowing()) {
		boolean updatingState = isUpdating;
		isUpdating = true;
		update_counter++;
		//DEBUG.printMessage(" *** Updating table: " + update_counter); // );
		//DEBUG.printStackTrace();

		// FIXME:gallery.updateSort();

		if (gallery.getSize() > 0)
		{
			// DEBUG.printMessage(String.format("updateTable: Table Model:\nsize: %d, cols: %d rows: %d",
			// gallery.getSize(), gallery.getColumnCount(),
			// gallery.getRowCount()));
			Dimension dimension = gallery.getCellSize();
			// renderer.setPreferredSize(dimension);

			// Adjusts rows size.
			table.setRowHeight(dimension.height);
			rowHeader.setFixedCellHeight(dimension.height);
			rowHeaderModel.setSize(gallery.getRowCount());
		}
		// Adjusts columns width
		if (gallery.adjustWidth)
			gallery.getColumnModel().adjustColumnsWidth(table);
		else
			gallery.adjustWidth = true;
		// columnModel.setWidth(dimension.width);

		// If auto adjust columns is enabled, refresh!
		jsRows.setValue(gallery.getRowCount());
		jsColumns.setValue(gallery.getColumnCount());

		rowHeader.revalidate();
		rowHeader.repaint();
		// table.revalidate();
		// repaint();
		// }

		jsZoom.setValue(data.zoom);
		isUpdating = updatingState;
		SwingUtilities.invokeLater(new Runnable()
		{
			public void run()
			{
				gallery.updateTableSelection(table);
			}
		});

		// }
	}// function updateTable

	/** Adjust the columns depending on the current windows width and cell width */
	private boolean adjustColumns()
	{
		if (autoAdjustColumns){
			//DEBUG.printStackTrace();
			return gallery.adjustColumn(width - 50);
		}
		return false;
		// DEBUG.printMessage(String.format(
		// "==>> JFrameGallery.autoAdjust: width: %d", width));
		// int rw = rowHeader.getWidth();
		// FIXME
		// gallery.autoAdjustColumns(
		// // jsPanel.getVisibleRect().width - rowHeader.getWidth(),
		// // jsPanel.getViewportBorderBounds().width -
		// // rowHeader.getWidth(),
		// // jspContent.getViewport().getWidth() - rowHeader.getWidth()
		// w - rw, table.getIntercellSpacing().width);
		// updateColumnsRowsValues();
	}

	public void setAutoAdjustColumns(boolean autoAdjustColumns)
	{
		this.autoAdjustColumns = autoAdjustColumns;
		jcbAutoAdjustColumns.setSelected(autoAdjustColumns);
		jsColumns.setEnabled(!autoAdjustColumns);
		jsRows.setEnabled(!autoAdjustColumns);
	}

	private void makeVisible(int index)
	{
		int coords[] = gallery.getCoords(index);
		//DEBUG.printMessage(String.format("gotoImage, index: %d, row: %d, col:%d", index, coords[0], coords[1]));

		// Gets current selected cell bounds.
		Rectangle rect = table.getCellRect(coords[0], coords[1], true);

		// Ensures item is visible
		Point pos = jspContent.getViewport().getViewPosition();
		rect.translate(-pos.x, -pos.y);
		jspContent.getViewport().scrollRectToVisible(rect);

		repaint();
	}

	private void goToImage(int index)
	{
		gallery.gotoItem(index);
		makeVisible(index);
	}

	public class Worker implements Runnable
	{
		public static final int STATS = 0;
		public static final int PCA = 1;
		public static final int FSC = 2;
		public String message;
		/** Constructor selecting operation */
		private int op; // store operation
		private MetaData imagesmd;

		public Worker(int operation, MetaData imagesmd)
		{
			op = operation;
			this.imagesmd = imagesmd;
			if (imagesmd.findObjects().length == 0)
				throw new IllegalArgumentException("No images available");
		}

		public void run()
		{
			try
			{
				switch (op)
				{
				case STATS:
					computeStatsImages(imagesmd);
					break;
				case PCA:
					pca(imagesmd);
					break;
				case FSC:
					fsc(imagesmd);
					break;
				}
				imagesmd.destroy();
			}
			catch (Exception e)
			{
				XmippWindowUtil.releaseGUI(GalleryJFrame.this.getRootPane());
				showException(e);
				return;

			}
			XmippWindowUtil.releaseGUI(GalleryJFrame.this.getRootPane());
		}

		public String getMessage()
		{
			switch (op)
			{
			case STATS:
				return "Computing average and std images...";
			case PCA:
				return "Computing PCA...";
			case FSC:
				return "Computing FSC...";
			}
			return "";
		}

	}

	/** Function to create and launch the worker, blocking the gui */
	public void runInBackground(int operation)
	{

		MetaData imagesmd = data.getImagesMd();
		Worker w = new Worker(operation, imagesmd);
		XmippWindowUtil.blockGUI(this, w.getMessage());
		Thread thr = new Thread(w);
		thr.start();

	}

	private void computeStatsImages(MetaData imagesmd) throws Exception
	{
		ImageGeneric imgAvg = new ImageGeneric();
		ImageGeneric imgStd = new ImageGeneric();

		imagesmd.getStatsImages(imgAvg, imgStd, data.useGeo, data.getRenderLabel());
		ImagePlus impAvg = XmippImageConverter.convertToImagePlus(imgAvg);
		ImagePlus impStd = XmippImageConverter.convertToImagePlus(imgStd);
		imgAvg.destroy();
		imgStd.destroy();

		XmippImageWindow winAvg = new XmippImageWindow(new ImagePlusLoader(impAvg), "AVG: " + data.getFileName());
		XmippWindowUtil.setLocation(0.2f, 0.5f, winAvg, this);
		winAvg.setVisible(true);
		XmippImageWindow winStd = new XmippImageWindow(new ImagePlusLoader(impStd), "STD: " + data.getFileName());

		XmippWindowUtil.setLocation(0.8f, 0.5f, winStd, this);
		winStd.setVisible(true);
		imagesmd.destroy();
	}

	private boolean openClassesDialog()
	{
		if (dlgClasses == null)
		{
			dlgClasses = new ClassesJDialog(GalleryJFrame.this);
		}
		boolean result = dlgClasses.showDialog();
		dlgClasses.resetClasses();
		return result;
	}

	static String forceExtension(String filename, String ext)
	{
		int dot = filename.lastIndexOf(".");
		return filename.substring(0, dot) + ext;
	}

	public void pca(MetaData imagesmd) throws Exception
	{
		ImageGeneric image = new ImageGeneric();
		imagesmd.getPCAbasis(image, data.getRenderLabel());
		ImagePlus imp = XmippImageConverter.convertToImagePlus(image);
		imp.setTitle("PCA: " + data.getFileName());
		ImagesWindowFactory.openXmippImageWindow(this, imp, false);
		imagesmd.destroy();

	}

	public void fsc(MetaData imagesmd) throws Exception
	{
		FSCJFrame frame = new FSCJFrame(data, imagesmd);
		XmippWindowUtil.centerWindows(frame, this);
		frame.setVisible(true);
	}

	/***
	 * Helper function to create toolbar toggle buttons
	 */
	protected void setupButton(AbstractButton btn, String icon, String text, ActionListener listener)
	{
		// Add toggle button to set/unset global normalization
		// JToggleButton btn = new javax.swing.JToggleButton();
		btn.setFocusable(false);
		btn.setIcon(XmippResource.getIcon(icon));
		btn.setToolTipText(text);
		btn.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
		btn.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
		btn.addActionListener(listener);
		// return btn;
	}

	private void updateViewState()
	{
		ImageIcon icon;
		String text;
		if (!data.isGalleryMode())
		{
			icon = XmippResource.VIEW_GALLERY_ICON;
			text = XmippLabel.LABEL_VIEW_GALLERY;
		}
		else
		{
			icon = XmippResource.VIEW_MD_ICON;
			text = XmippLabel.LABEL_VIEW_MD;
		}
		btnChangeView.setIcon(icon);
		btnChangeView.setToolTipText(text);
		boolean allowColsResize = true;
		if (data.isTableMode())
		{ // if we are in table mode only allow change
			// if exist render label
			boolean hasRender = data.allowGallery();
			btnChangeView.setEnabled(hasRender);
			jsZoom.setEnabled(hasRender);
			jlZoom.setEnabled(hasRender);
			boolean isCol = data.isColumnFormat();
			allowColsResize = false;
			jsGoToImage.setEnabled(isCol && gallery.getSize() > 0);
			jlGoToImage.setEnabled(isCol);
		}
		jsColumns.setEnabled(allowColsResize);
		jlColumns.setEnabled(allowColsResize);
		jsRows.setEnabled(allowColsResize);
		jlRows.setEnabled(allowColsResize);
		jcbAutoAdjustColumns.setEnabled(allowColsResize);
	}

	public void reloadTableData()
	{
		reloadTableData(true);
	}

	/** Reload table data */
	public void reloadTableData(boolean changed)
	{
		try
		{
			//DEBUG.printMessage("reloadTableData...");
			if (table != null)
				table.removeAll();
			createModel();
			// gallery.setShowLabels(menu.getShowLabel());
			createTable();

			menu.update();
			updateVisibleCombos();
			if (dlgSave != null && changed)
				dlgSave.setInitialValues();

			this.saved = !changed;

			setGalleryTitle();

			if (rows != null)
				gallery.setRows(rows);
			if (columns != null)
				gallery.setColumns(columns);

		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	private void reloadMd() throws Exception
	{
		reloadMd(true);
	}

	/**
	 * Reload metadata info, rebuild the table This function is called whenever
	 * a change is made on metadata, that's why changes are reported to
	 * GalleryData
	 * */

	private void reloadMd(boolean changed) throws Exception
	{
		data.loadMd();
		reloadTableData(changed);
		data.setMdChanges(changed);

	}// function reloadMd

	/**
	 * Fill some label mode can be: "constant", "linear", "uniform", "gaussian"
	 * values is a list of string
	 * */
	public void fillLabel(int label, String mode, String... values) throws Exception
	{
		if (mode.equalsIgnoreCase(MetaData.FILL_CONSTANT))
			data.md.fillConstant(label, values[0]);
		else
		{
			Double v1 = Double.parseDouble(values[0]);
			Double v2 = Double.parseDouble(values[1]);
			if (mode.equalsIgnoreCase(MetaData.FILL_LINEAR))
				data.md.fillLinear(label, v1, v2);
			else if (mode.equalsIgnoreCase(MetaData.FILL_RAND_UNIFORM))
				data.md.fillRandom(label, "uniform", v1, v2);
			else if (mode.equalsIgnoreCase(MetaData.FILL_RAND_GAUSSIAN))
				data.md.fillRandom(label, "gaussian", v1, v2);
			else
				throw new Exception("Unknown label fill mode: " + mode);
		}
		reloadMd();
	}

	/**
	 * Delete selected or disabled items if 'selected' is true, selection is
	 * removed if false, the disabled items
	 * */
	public void removeObjects(boolean selected)
	{
		try
		{
			String type = selected ? "selected" : "disabled";
			if (XmippDialog.showWarning(this, String.format("Are you sure to delete %s items?", type)))
			{
				if (selected)
					data.removeSelection();
				else
					data.md.removeDisabled();
				reloadMd();
			}
		}
		catch (Exception ex)
		{
			showException(ex);
		}
	}

	/** Save selected items as a metadata */
	public void saveSelection() throws Exception
	{
		MetaData md = data.getSelectionMd();
		SaveJDialog dlg = new SaveJDialog(this, "selection.xmd", true);
		boolean save = dlg.showDialog();
		if (save)
		{
			boolean overwrite= dlg.isOverwrite();
			String path = dlg.getMdFilename();
			saveSelection(path, overwrite);
		}
		md.destroy();
	}
        
        	/** Save selected items as a metadata */
	public void saveSelection(String path, boolean overwrite) throws Exception
	{
		MetaData md = data.getSelectionMd();
                String file = path.substring(path.lastIndexOf("@") + 1, path.length());
                if (!new File(file).exists())// overwrite or append, save selection
                        md.write(path);
                else
                {
                        if (overwrite)
                                md.write(path);// overwrite with active block only, other
                                                                // blocks were dismissed
                        else
                                md.writeBlock(path);// append selection

                }
		
		md.destroy();
	}

	/** Find and replace in metadata */
	public void findReplace() throws Exception
	{
		MDSearchJDialog dlg = new MDSearchJDialog(this, table, data.md);
		dlg.setVisible(true);
	}

	/** Drop some label from the metadata */
	public void removeLabel(int label) throws Exception
	{
		data.md.removeLabel(label);
		reloadMd();
	}

	/***
	 * Function to create the main toolbar
	 */
	protected void createToolbar()
	{
		// Create Main TOOLBAR
		toolBar = new JToolBar();
		toolBar.setRollover(true);
		toolBar.setLayout(new FlowLayout(FlowLayout.LEFT));

		btnChangeView = new JButton();
		// updateViewState();
		btnChangeView.addActionListener(new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				data.changeMode();
				reloadTableData();
				makeVisible(gallery.getFirstSelectedIndex());
			}
		});

		toolBar.add(btnChangeView);
		toolBar.addSeparator();

		jlZoom = new javax.swing.JLabel();
		jsZoom = new javax.swing.JSpinner();
		jlZoom.setIcon(XmippResource.getIcon(XmippResource.ZOOM));
		jlZoom.setToolTipText(XmippLabel.LABEL_ZOOM);
		toolBar.add(jlZoom);

		jsZoom.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(1), Integer.valueOf(1), null, Integer.valueOf(1)));
		jsZoom.addChangeListener(new javax.swing.event.ChangeListener()
		{
			public void stateChanged(javax.swing.event.ChangeEvent evt)
			{
				Integer zoom = (Integer) jsZoom.getValue();
				if (gallery.getCellSize().getHeight() < 30)
				{
					jsZoom.setValue(gallery.data.zoom);//keep previous zoom
					return;
				}
				gallery.setZoom(zoom);
				makeVisible(gallery.getFirstSelectedIndex());
				// gallery.updateTableSelection(table);
			}
		});

		toolBar.add(jsZoom);
		toolBar.addSeparator();

		jlGoToImage = new javax.swing.JLabel();
		jsGoToImage = new javax.swing.JSpinner();
		jlGoToImage.setIcon(XmippResource.getIcon(XmippResource.GOTO));
		jlGoToImage.setToolTipText(XmippLabel.LABEL_GOTO_ITEM);
		toolBar.add(jlGoToImage);

		if (gallery.getSize() > 0)
			jsGoToImage.setValue(1);

		jsGoToImage.addChangeListener(new javax.swing.event.ChangeListener()
		{
			public void stateChanged(javax.swing.event.ChangeEvent evt)
			{
				jsGoToImageStateChanged(evt);
			}
		});
		toolBar.add(jsGoToImage);

		toolBar.addSeparator();

		jcbAutoAdjustColumns = new JToggleButton();
		setupButton(jcbAutoAdjustColumns, XmippResource.ADJUST_COLS, XmippLabel.MSG_ADJUST_COLS, new java.awt.event.ActionListener()
		{
			public void actionPerformed(java.awt.event.ActionEvent evt)
			{
				autoAdjustColumns(jcbAutoAdjustColumns.isSelected());
			}
		});
		jcbAutoAdjustColumns.setSelected(true);
		jlRows = new javax.swing.JLabel();
		jsRows = new javax.swing.JSpinner();
		jlColumns = new javax.swing.JLabel();
		jsColumns = new javax.swing.JSpinner();
		toolBar.add(jcbAutoAdjustColumns);

		jlColumns.setText(XmippLabel.LABEL_COLUMNS);
		toolBar.add(jlColumns);

		jsColumns.addChangeListener(new javax.swing.event.ChangeListener()
		{
			public void stateChanged(javax.swing.event.ChangeEvent evt)
			{
				jsColumnsStateChanged(evt);
			}
		});
		toolBar.add(jsColumns);

		jlRows.setText(XmippLabel.LABEL_ROWS);
		toolBar.add(jlRows);

		jsRows.addChangeListener(new javax.swing.event.ChangeListener()
		{
			public void stateChanged(javax.swing.event.ChangeEvent evt)
			{
				jsRowsStateChanged(evt);
			}
		});
		toolBar.add(jsRows);

		// Some settings of the spinners

		if (gallery.getSize() > 0)
		{

			jsRows.setModel(new SpinnerNumberModel(1, 1, gallery.getSize(), 1));
			jsColumns.setModel(new SpinnerNumberModel(1, 1, gallery.getSize(), 1));
			jsGoToImage.setModel(new SpinnerNumberModel(1, 1, gallery.getSize(), 1));
		}

		int TEXTWIDTH = 4;
		((JSpinner.NumberEditor) jsZoom.getEditor()).getTextField().setColumns(TEXTWIDTH);
		((JSpinner.NumberEditor) jsGoToImage.getEditor()).getTextField().setColumns(TEXTWIDTH);
		((JSpinner.NumberEditor) jsRows.getEditor()).getTextField().setColumns(TEXTWIDTH);
		((JSpinner.NumberEditor) jsColumns.getEditor()).getTextField().setColumns(TEXTWIDTH);
		
		initResliceButtonMenu();
		toolBar.add(reslicebt);

	}// function createToolbar

	/** Create combos for selection of block and volume if its the case */
	protected void createCombos()
	{
		cbPanel = new JPanel();
		cbPanel.setLayout(new FlowLayout(FlowLayout.LEFT));

		// Add blocks selector combo
		jlBlocks = new JLabel(XmippLabel.LABEL_BLOCK);
		jcbBlocks = new JComboBox();
		DEBUG.printFormat("data.getNumberOfBlocks: %d\n", data.getNumberOfBlocks());
		
		if (data.getNumberOfBlocks() > 0)
		{
			cbPanel.add(jlBlocks);
			jcbBlocks.setModel(new ComboBoxModel()
			{

				@Override
				public int getSize()
				{
					return data.getNumberOfBlocks();
				}

				@Override
				public Object getElementAt(int index)
				{
					return data.getBlock(index);
				}

				@Override
				public void setSelectedItem(Object item)
				{
					if (proceedWithChanges())
					{
						data.selectBlock((String) item);
						jcbVolumes.invalidate();
						try
						{
							data.loadMd();
							reloadTableData();
						}
						catch (Exception e)
						{
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}

				@Override
				public Object getSelectedItem()
				{
					return data.selectedBlock;
				}

				@Override
				public void removeListDataListener(ListDataListener arg0)
				{
				}

				@Override
				public void addListDataListener(ListDataListener arg0)
				{
					// TODO Auto-generated method stub
				}
			});

			cbPanel.add(jcbBlocks);
		}
		// Add volumes selector combo
		jlVolumes = new JLabel(XmippLabel.LABEL_VOLUME);
		cbPanel.add(jlVolumes);
		jcbVolumes = new JComboBox();
		jcbVolumes.setModel(new ComboBoxModel()
		{

			@Override
			public int getSize()
			{
				return data.getNumberOfVols();
			}

			@Override
			public Object getElementAt(int index)
			{
				return removePrefix(data.getVolumeAt(index));
			}

			@Override
			public void setSelectedItem(Object anItem)
			{
				data.selectVolume(data.commonVolPrefix + (String) anItem);
				reloadTableData();
			}

			@Override
			public Object getSelectedItem()
			{
				return removePrefix(data.selectedVolFn);
			}

			@Override
			public void addListDataListener(ListDataListener arg0)
			{
				// TODO Auto-generated method stub

			}

			@Override
			public void removeListDataListener(ListDataListener arg0)
			{
				// TODO Auto-generated method stub

			}

			public String removePrefix(String value)
			{
				return value.replaceFirst(data.commonVolPrefix, "");
			}
		});
		cbPanel.add(jcbVolumes);
	}

	protected void updateVisibleCombos()
	{

		boolean showBlocks = data.getNumberOfBlocks() > 0;
		boolean showVols = data.getNumberOfVols() > 1 && data.isVolumeMode();
		jcbBlocks.setVisible(showBlocks);
		jcbVolumes.setVisible(showVols);
		jlBlocks.setVisible(showBlocks);
		jlVolumes.setVisible(showVols);
		cbPanel.setVisible(showBlocks || showVols);
		//pack();
	}

	private void jsRowsStateChanged(javax.swing.event.ChangeEvent evt)
	{// GEN-FIRST:event_jsRowsStateChanged
		if (!isUpdating)
		{
			rows = (Integer) jsRows.getValue();
			columns = null;
			gallery.setRows(rows);
		}
	}

	private void jsColumnsStateChanged(javax.swing.event.ChangeEvent evt)
	{// GEN-FIRST:event_jsColumnsStateChanged
		if (!isUpdating)
		{
			columns = (Integer) jsColumns.getValue();
			rows = null;
			gallery.setColumns(getLastColums());

		}
	}

	public int getLastRows()
	{
		return rows;
	}

	public int getLastColums()
	{
		return columns;
	}

	private void jsGoToImageStateChanged(javax.swing.event.ChangeEvent evt)
	{
		if (!isUpdating)
		{
			Integer intValue = (Integer) jsGoToImage.getValue();
			if (intValue < 0)
				intValue = 1;
			else if (intValue >= gallery.getSize())
				intValue = gallery.getSize();
			goToImage(intValue - 1);
		}
	}

	private void formComponentResized(java.awt.event.ComponentEvent evt)
	{
		width = getSize().width;
		if (!isUpdating && autoAdjustColumns)
		{
			adjustColumns();
		}
	}

	public void selectItem(int row, int col)
	{
		if (row < 0 || row > table.getRowCount() - 1 || col < 0 || col > table.getColumnCount() - 1)
			return;
		
		if (gallery.data.isGalleryMode() && 
			row * table.getColumnCount() + col + 1 > gallery.getSize())
		{
			int[] coords = gallery.getCoords(gallery.getSize() - 1);
			row = coords[0];
			col = coords[1];

		}

		gallery.clearSelection();
		gallery.touchItem(row, col);
		makeVisible(row);

	}

	protected void tableMouseClicked(MouseEvent evt)
	{
		final Point p = evt.getPoint();
		int row = table.rowAtPoint(p);
		int col = table.columnAtPoint(p);

		if (evt.getButton() == MouseEvent.BUTTON1)
		{ // Left click.
			if (evt.getClickCount() > 1)
			{
				try
				{
					gallery.handleDoubleClick(row, col);
				}
				catch (Exception e)
				{
					XmippDialog.showError(this, e.getMessage());
				}
			}
			else
			{
				// Ctrl adds items to selection, otherwise previous ones are
				// removed.
				boolean move = true;
				if (!evt.isControlDown() && !evt.isShiftDown())
				{

					boolean clear = true;

					if (gallery.getSelectionCount() > 1)
					{
						clear = XmippDialog.showQuestion(this, "You will lose previous selection.\nDo you want to proceed?");
					}

					if (clear)
					{
						gallery.clearSelection();
						gallery.touchItem(row, col);

					}
					else
					{
						gallery.fireTableDataChanged();
						move = false;
					}
				}
				else
				{
					if (evt.isShiftDown())
					{
						gallery.selectRange(previousSelectedRow, previousSelectedCol, row, col, true);
					}
					else if (evt.isControlDown())
					{
						gallery.touchItem(row, col);
					}
				}
				if (move)
				{
					isUpdating = true;
					int index = gallery.getIndex(row, col);
					if (gallery.isValidIndex(index))
						jsGoToImage.setValue(index + 1);
					isUpdating = false;
				}
				if (!evt.isShiftDown())
				{
					previousSelectedRow = row;
					previousSelectedCol = col;
				}
			}

		}
		else if (evt.getButton() == MouseEvent.BUTTON3)
		{ // Right click.

			final MouseEvent me = evt;
			if (gallery.handleRightClick(row, col, jpopUpMenuTable))
			{
				if (gallery.getSelectionCount() < 2)
				{
					gallery.clearSelection();
					gallery.touchItem(row, col);
				}
				SwingUtilities.invokeLater(new Runnable()
				{
					public void run()
					{
						jpopUpMenuTable.show(me.getComponent(), p);
					}
				});
			}
		}
		table.invalidate();
		table.repaint();
		refreshExtractFrame();
	}// function tableMouseClicked

	private void autoAdjustColumns(boolean value)
	{
		setAutoAdjustColumns(value);
		adjustColumns();
	}

	private void setResliceView(int view)
	{
		data.resliceView = view;
		reloadTableData();
	}

	class GalleryMenu extends XmippMenuBarCreator
	{

		private QuickHelpJDialog quickhelpdlg;

		@Override
		protected void createItems() throws Exception
		{
			// File
			addItem(FILE, "File");
			addItem(FILE_OPEN, "Open ...", null, "control released O");
			addItem(FILE_OPENWITH_IJ, "Open with ImageJ", "ij.gif", "control released J");
			addItem(FILE_OPENWITH_CHIMERA, "Open with Chimera", "chimera.gif", "control released H");
			addItem(FILE_OPENMICROGRAPHS, "Open Particle Micrographs");
			addItem(FILE_INFO, "File info ...");

			addSeparator(FILE);
			addItem(FILE_SAVE, "Save", "save.gif", "control released S");
			addItem(FILE_SAVEAS, "Save as", "save_as.gif");
			addItem(FILE_REFRESH, "Refresh", "refresh.gif", "released F5");
			addSeparator(FILE);
			addItem(FILE_EXIT, "Exit", null, "control released Q");
			// Display
			addItem(DISPLAY, "Display");
			addItem(DISPLAY_NORMALIZE, "Global normalization", null, "control released N");
			addItem(DISPLAY_SHOWLABELS, "Display labels", null, "control released L");
			addSeparator(DISPLAY);
			addItem(DISPLAY_RENDERIMAGES, "Render images", null, "control released R");
			
			addRenderImageColumnItem();
			
		
			addItem(DISPLAY_APPLYGEO, "Apply geometry", null, "control released G");
			addItem(DISPLAY_WRAP, "Wrap", null, "control released W");
			addItem(DISPLAY_COLUMNS, "Columns ...", "columns.gif");
			addItem(DISPLAY_RESLICE, "Reslice");
			for (int i = 0; i < ImageGeneric.VIEWS.length; ++i)
				addItem(DISPLAY_RESLICE_VIEWS[i], reslices[i]);
			// Metadata operations
			addItem(METADATA, "Metadata");
			addItem(STATS, "Statistics");
			addItem(STATS_AVGSTD, "Avg & Std images");
			addItem(STATS_PCA, "PCA");
			addItem(STATS_FSC, "FSC");
			addItem(MD_PLOT, "Plot", "plot.png");
			addItem(MD_CLASSES, "Classes");
			addItem(MD_EDIT_COLS, "Edit labels", "edit.gif");
			addItem(MD_ADD_OBJECT, "Add new object", "new_object.gif");
			addItem(MD_REMOVE_DISABLED, "Remove disabled", "delete.gif");
			addItem(MD_REMOVE_SELECTION, "Remove selection");
			addItem(MD_SAVE_SELECTION, "Save selection", "save.gif");
			addSeparator(METADATA);
			addItem(MD_FIND_REPLACE, "Find & Replace", "search.gif", "control released F");
			// Help
			addItem(HELP, "Help");
			addItem(HELP_ONLINE, "Online help", "online_help.gif");
			addItem(KEY_ASSIST, "Key Assist...");
		}// function createItems

		public void update()
		{
			boolean galMode = data.isGalleryMode();
			boolean volMode = data.isVolumeMode();
			setItemEnabled(FILE_OPENWITH_CHIMERA, volMode);
			setItemEnabled(FILE_OPENMICROGRAPHS, data.hasMicrographParticles());
			setItemEnabled(FILE_SAVE, !volMode);
			setItemEnabled(FILE_SAVEAS, !volMode);
			setItemSelected(DISPLAY_NORMALIZE, gallery.getNormalized());
			setItemEnabled(DISPLAY_APPLYGEO, data.containsGeometryInfo());
			setItemEnabled(DISPLAY_WRAP, data.containsGeometryInfo() && data.useGeo);
			setItemSelected(DISPLAY_WRAP, data.wrap);
			setItemSelected(DISPLAY_APPLYGEO, data.useGeo);
			setItemEnabled(DISPLAY_RENDERIMAGES, !galMode && data.hasRenderLabel());
			setItemSelected(DISPLAY_RENDERIMAGES, data.globalRender);
			setItemEnabled(DISPLAY_RENDERIMAGECOLUMN, galMode);
			for (int i = 0; i < ImageGeneric.VIEWS.length; ++i)
				setItemSelected(DISPLAY_RESLICE_VIEWS[i], (data.resliceView == ImageGeneric.VIEWS[i]));
			setItemEnabled(DISPLAY_COLUMNS, !galMode);
			setItemEnabled(DISPLAY_RESLICE, volMode);
			setItemEnabled(MD_CLASSES, data.is2DClassificationMd());
			setItemEnabled(MD_PLOT, data.isTableMode());
			boolean isCol = data.isColumnFormat();
			setItemEnabled(STATS, isCol && !volMode);
			setItemEnabled(MD_ADD_OBJECT, isCol);
			setItemEnabled(MD_REMOVE_DISABLED, isCol);
			setItemEnabled(MD_REMOVE_SELECTION, isCol);
			setItemEnabled(MD_SAVE_SELECTION, isCol);
			setItemEnabled(MD_FIND_REPLACE, isCol && !galMode);
			reslicebt.setEnabled(volMode);
		}// function update

		@Override
		protected void handleActionPerformed(ActionEvent evt)
		{
			String cmd = evt.getActionCommand();
			try
			{
				if (cmd.equals(DISPLAY_NORMALIZE))
				{
					gallery.setNormalized(getItemSelected(DISPLAY_NORMALIZE));
				}
				else if (cmd.equals(DISPLAY_APPLYGEO) || cmd.equals(DISPLAY_WRAP))
				{
					if (data.containsGeometryInfo())
					{
						((MetadataGalleryTableModel) gallery).setUseGeometry(getItemSelected(DISPLAY_APPLYGEO), getItemSelected(DISPLAY_WRAP));
						setItemEnabled(DISPLAY_WRAP, data.containsGeometryInfo() && data.useGeo);

					}
				}
				else if (cmd.equals(DISPLAY_SHOWLABELS))
				{
					gallery.setShowLabels(getItemSelected(DISPLAY_SHOWLABELS));
				}
				else if (cmd.equals(DISPLAY_RENDERIMAGES))
				{
					gallery.setRenderImages(getItemSelected(DISPLAY_RENDERIMAGES));
					makeVisible(gallery.getFirstSelectedIndex());
				}
				
				else if (cmd.equals(DISPLAY_COLUMNS))
				{
					ColumnsJDialog dialog = new ColumnsJDialog(GalleryJFrame.this);
					boolean result = dialog.showDialog();
					if (result)
					{
						ArrayList<ColumnInfo> columns = dialog.getColumnsResult();
						isUpdating = true;
						((MetadataGalleryTableModel) gallery).updateColumnInfo(columns);
						gallery.fireTableDataChanged();
						setItemEnabled(DISPLAY_RENDERIMAGES, data.globalRender);
						// menu.enableRenderImages(data.globalRender);
						isUpdating = false;
					}
				}
				else if (cmd.equals(STATS_AVGSTD))
					runInBackground(Worker.STATS);
				else if (cmd.equals(STATS_PCA))
					runInBackground(Worker.PCA);
				else if (cmd.equals(STATS_FSC))
					runInBackground(Worker.FSC);
				else if (cmd.equals(FILE_OPEN))
				{
					if (fc.showOpenDialog(GalleryJFrame.this) != XmippFileChooser.CANCEL_OPTION)
					{
						if (Filename.exists(fc.getSelectedFile().getPath()))
							ImagesWindowFactory.openFileAsDefault(fc.getSelectedPath());
						else
							XmippDialog.showError(GalleryJFrame.this, String.format("File: '%s' doesn't exist.", fc.getSelectedPath()));
					}
				}
				else if (cmd.equals(FILE_SAVE))
				{
					save();
				}
				else if (cmd.equals(FILE_SAVEAS))
				{
					saveAs();
				}
				else if (cmd.equals(FILE_EXIT))
				{
					close();
				}
				else if (cmd.equals(FILE_OPENWITH_CHIMERA))
				{
					try
					{
						String args = data.selectedVolFn;
						// if (Filename.isSpiderVolume(args))
						// args = "spider:" + args;
						// // FIXME: Check chimera is installed
						// Process p = new ProcessBuilder("chimera",
						// args).start();
						Runtime.getRuntime().exec("xmipp_chimera_client -i " + args);
					}
					catch (Exception ex)
					{
						ex.printStackTrace();
					}
				}
				else if (cmd.equals(FILE_OPENMICROGRAPHS))
				{
					openMicrographs();
				}

				else if (cmd.equals(FILE_INFO))
				{
					XmippDialog.showInfo(GalleryJFrame.this, data.getFileInfo());
				}

				else if (cmd.equals(FILE_OPENWITH_IJ))
				{
					try
					{
						ImagePlusLoader loader = gallery.getImageLoader();
						ImagesWindowFactory.openXmippImageWindow(GalleryJFrame.this, loader, true);
						XmippIJUtil.showImageJ(Tool.VIEWER);
					}
					catch (Exception e1)
					{
						e1.printStackTrace();
					}
				}
				else if (cmd.equals(FILE_REFRESH))
				{
					data.readMd();
					reloadTableData();

				}
				else if (cmd.contains(DISPLAY_RESLICE))
				{
					for (int i = 0; i < ImageGeneric.VIEWS.length; ++i)
						if (cmd.equals(DISPLAY_RESLICE_VIEWS[i]))
						{
							setResliceView(ImageGeneric.VIEWS[i]);
							break;
						}
				}
				else if (cmd.equals(MD_PLOT))
				{
					PlotJDialog dlg = new PlotJDialog(GalleryJFrame.this);
					dlg.showDialog();
				}
				else if (cmd.equals(MD_CLASSES))
				{
					openClassesDialog();
				}
				else if (cmd.equals(MD_EDIT_COLS))
				{
					EditLabelsJDialog dlg = new EditLabelsJDialog(GalleryJFrame.this);
					dlg.showDialog();
				}
				else if (cmd.equals(MD_REMOVE_SELECTION))
				{
					if (gallery.getSelectionCount() > 0)
						removeObjects(true);
				}
				else if (cmd.equals(MD_REMOVE_DISABLED))
				{
					removeObjects(false);
				}
				else if (cmd.equals(MD_SAVE_SELECTION))
				{
					saveSelection();
				}
				else if (cmd.equals(MD_FIND_REPLACE))
				{
					findReplace();
				}
				else if (cmd.equals(MD_ADD_OBJECT))
				{
					AddObjectJDialog dlg = new AddObjectJDialog(GalleryJFrame.this);
					if (dlg.showDialog())
					{
						data.md.unionAll(dlg.md);
						reloadMd();
					}
				}
				else if (cmd.equals(HELP_ONLINE))
				{
					XmippWindowUtil.openURI("http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/WebHome");
				}
				else if (cmd.equals(KEY_ASSIST))
				{
					if (quickhelpdlg == null)
						quickhelpdlg = new QuickHelpJDialog(GalleryJFrame.this, false, "Key Assist", getKeyAssist());
					quickhelpdlg.setVisible(true);

				}

			}
			catch (Exception e)
			{
				showException(e);
			}
		}// function handleActionPerformed
		
		
		protected void addRenderImageColumnItem()
		{                  
			addItem(DISPLAY_RENDERIMAGECOLUMN, "Render Image Column");
			JMenuItem mi;
			
			imagecolumns = new Hashtable<String, ColumnInfo>();
			for (ColumnInfo column : data.getColumns())
				if (column.allowRender)
					imagecolumns.put(column.toString(), column);
			boolean rendercolumn = imagecolumns.size() > 1;
			setItemEnabled(DISPLAY_RENDERIMAGECOLUMN, rendercolumn);
			if(rendercolumn)
			{
				// Create the popup menu.
				String id, column;
				Enumeration<String> keys = imagecolumns.keys();
				while (keys.hasMoreElements())
				{
					column = keys.nextElement();
					id = String.format("Display.RenderImagesColumn.%s_rb", column);
					mi = addItem(id, column);
					mi.addActionListener(new RenderColumnActionListener());
					if(data.getRenderColumn().toString().equals(column))
						setItemSelected(id, true);
				}
			}

		}
		
		class RenderColumnActionListener implements ActionListener
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				JRadioButtonMenuItem mi = (JRadioButtonMenuItem)e.getSource();
				String key = mi.getText();
				data.setRenderColumn(imagecolumns.get(key));
				reloadTableData();
			}

		}
			
	}// class GalleryMenu //////////////////////////////////////////////////////////
	
	
	
	
	
	
	
	
	

	class GalleryPopupMenu extends XmippPopupMenuCreator
	{
		protected int row;
		protected int col;

		@Override
		protected void createItems() throws Exception
		{
			addItem(ENABLED, "Enable", "enable.gif");
			addItem(DISABLED, "Disable", "disable.gif");
			addItem(REFRESH, "Refresh", "refresh.gif");
			// addSeparator();
			addItem(OPEN, "Open");
			addItem(OPEN_ASTEXT, "Open as text");
			addItem(CTF_PROFILE, "Show CTF profile");
			addItem(CTF_RECALCULATE, "Recalculate CTF");
			addSeparator();
			addItem(OPEN_IMAGES, "Open images");
			addItem(SAVE_IMAGES, "Save images", "save.gif");
			addItem(SET_CLASS, "Set class");
			addItem(SELECT, "Select");
			addItem(SELECT_ALL, "All", null, "control released A");
			addItem(SELECT_TOHERE, "To here");
			addItem(SELECT_FROMHERE, "From here");
			initItems();
		}// function createItems

		public void show(Component cmpnt, Point location)
		{
			setItemVisible(SET_CLASS, data.is2DClassificationMd());
			// This item visibility depends on current selection
			setItemVisible(SAVE_IMAGES, data.is2DClassificationMd() && gallery.getSelectionCount() > 0);
			setItemVisible(OPEN_IMAGES, data.is2DClassificationMd() && gallery.getSelectionCount() == 1);
					
			// Update menu items status depending on item.
			row = table.rowAtPoint(location);
			col = table.columnAtPoint(location);
			getPopupMenu().show(cmpnt, location.x, location.y);

		}// function show

		private void selectRange(int first, int last)
		{
			gallery.selectRange(first, last, true);
		}

		private void showCTF(boolean profile)
		{
			try
			{
				String ctfModel = data.md.getValueString(MDLabel.MDL_CTF_MODEL, data.ids[row]);
				String displayFilename = data.md.getValueString(MDLabel.MDL_PSD_ENHANCED, data.ids[row]);
				String psdFile = data.md.getValueString(MDLabel.MDL_PSD, data.ids[row]);

				ImageGeneric img = new ImageGeneric(displayFilename);
				ImagePlus imp = XmippImageConverter.readToImagePlus(img);

				if (profile)
					ImagesWindowFactory.openCTFWindow(imp, ctfModel, psdFile);
				else
				{
					MetaData mdRow = new MetaData(); 
					MDRow row2 = new MDRow();
					data.md.getRow(row2, data.ids[row]);
					mdRow.setRow(row2, mdRow.addObject());
					String sortFn = psdFile.replace(".psd", ".xmd");
					mdRow.write(sortFn);
					mdRow.destroy();
					ImagesWindowFactory.openCTFImage(imp, ctfModel, psdFile, ctfTasks, data.md.getFilename(), row, sortFn);
				}

			}
			catch (Exception e)
			{
				XmippDialog.showError(GalleryJFrame.this, e.getMessage());
			}
		}

		/** Set values to defaults */
		@Override
		public void initItems()
		{
			setItemVisible(OPEN, false);
			setItemVisible(OPEN_ASTEXT, false);
			setItemVisible(CTF_PROFILE, false);
			setItemVisible(CTF_RECALCULATE, false);
		}

		@Override
		protected void handleActionPerformed(ActionEvent evt)
		{
			String cmd = evt.getActionCommand();
			if (cmd.equals(SELECT_ALL))
			{
				selectRange(0, gallery.getSize() - 1);
			}
			else if (cmd.equals(SELECT_TOHERE))
			{
				selectRange(0, gallery.getIndex(row, col));
			}
			else if (cmd.equals(SELECT_FROMHERE))
			{
				selectRange(gallery.getIndex(row, col), gallery.getSize() - 1);
			}
			else if (cmd.equals(ENABLED))
			{
				gallery.setSelectionEnabled(true);
				// gallery.clearSelection();
				refreshExtractFrame();

			}
			else if (cmd.equals(DISABLED))
			{
				gallery.setSelectionEnabled(false);
				// gallery.clearSelection();
				refreshExtractFrame();

			}
			else if (cmd.equals(REFRESH))
			{
				gallery.refreshAt(row, col);

			}
			else if (cmd.equals(OPEN))
			{
                               
                                MetadataGalleryTableModel mg = (MetadataGalleryTableModel) gallery;

                                ColumnInfo ci = mg.visibleLabels.get(col);

                                if (ci.allowRender)
                                        gallery.handleDoubleClick(row, col);
                                else
                                {
                                        String file = data.getValueFromCol(row, ci);

                                        ImagesWindowFactory.openFileAsDefault(file);
                                }
                                
			}
			else if (cmd.equals(OPEN_ASTEXT))
			{
				String file = gallery.getValueAt(row, col).toString();
				ImagesWindowFactory.openFileAsText(file, null);
			}
			else if (cmd.equals(CTF_PROFILE))
			{
				showCTF(true);
			}
			else if (cmd.equals(CTF_RECALCULATE))
			{
				showCTF(false);
			}
			else if (cmd.equals(SET_CLASS))
			{
				if (openClassesDialog())
				{
					int classNumber = dlgClasses.getSelectedClass();
					// DEBUG.printMessage(String.format("class: %d",
					// classNumber));
					gallery.setSelectionClass(classNumber);
				}
			}
			else if (cmd.equals(OPEN_IMAGES))
			{
				int index = gallery.getIndex(row, col);
				MetaData md = data.getClassImages(index);
				if (md != null)
					openMetadata(md);
				else
					XmippDialog.showWarning(GalleryJFrame.this, "This class has no images");
			}
			else if (cmd.equals(SAVE_IMAGES))
			{
				File f = new File(data.getFileName());
				SaveImagesJDialog dialog = new SaveImagesJDialog(GalleryJFrame.this, f.getParent() + "/images_selection.xmd");
				dialog.showDialog();					
			}
			initItems();

		}

	}// class JPopUpMenuGallery

	public void showException(Exception e)
	{
		XmippDialog.showException(this, e);
	}

	@Override
	public void setRunning(boolean running)
	{
		// XmippDialog.showInfo(this, String.format("Calculating ctf"));
		// XmippWindowUtil.blockGUI(getRootPane(), "Calculating CTF");
	}

	@Override
	public void setRowBusy(int row)
	{

		((MicrographsTableModel) gallery).setRowBusy(row);
	}

	@Override
	public void setRowIdle(int row)
	{
		((MicrographsTableModel) gallery).setRowIdle(row);

	}

	

	@Override
	public void done()
	{
		XmippDialog.showInfo(this, String.format("Calculating ctf: DONE"));
	}

	private void saveMd() throws Exception
	{
		saveMd(dlgSave.getMdFilename(), false);
	}

	private void saveMd(String path, boolean saveall) throws Exception
	{
		try
		{
			if (path == null)
				throw new IllegalArgumentException();

			boolean overwritewithblock;
			String file;
			if (path.contains("@"))
				file = path.substring(path.lastIndexOf("@") + 1, path.length());
			else
			{
				file = path;
				path = getBlock() + "@" + file;
			}

			File iofile = new File(file);
			if (!iofile.exists())// overwrite or append, save active
			{
				if (iofile.getParentFile() != null)
					iofile.getParentFile().mkdirs();
				data.md.write(path);
			}
			else
			{
				overwritewithblock = dlgSave.isOverwrite() && dlgSave.saveActiveMetadataOnly();
				if (overwritewithblock)
					data.md.write(path);// overwrite with active block only,
										// other blocks were dismissed
				else
					data.md.writeBlock(path);// either if save active block or all, save active, other blocks where already managed

			}
			if (!saveall)
			{
				data.setMdChanges(false);
				gallery.data.setFileName(file);
				if (path.contains("@"))
					gallery.data.selectBlock(path.substring(0, path.lastIndexOf("@")));
				reloadFile(file, false);
				setTitle(gallery.getTitle());
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}// function saveMd

        
	public String getFilename()
	{
		return data.getFileName();
	}
        
        
	private void saveAll() throws Exception
	{
		String from = data.getFileName();
		String blockto = dlgSave.getMdFilename();
		String to;
		if (blockto.contains("@"))
			to = blockto.substring(blockto.lastIndexOf("@") + 1, blockto.length());
		else
		{
			to = blockto;
			blockto = getBlock() + "@" + blockto;
		}
		
		if (from != null)
		{
			MetaData md;
			Hashtable<String, MetaData> mds = new Hashtable<String, MetaData>();
			for (String blockit : data.mdBlocks)
				mds.put(blockit, new MetaData(blockit + "@" + from));
			File file = new File(to);
			if (dlgSave.isOverwrite())
				file.delete();
			if (!file.exists() && file.getParentFile() != null)
				file.getParentFile().mkdirs();
			for (String blockit : data.mdBlocks)
			{
				md = mds.get(blockit);
				if (blockit.equals(getBlock()))
					saveMd(blockto, true);
				else
					md.writeBlock(blockit + "@" + to);
				md.destroy();
			}
		}
		else {
			saveMd(blockto, true);
		}

		data.setMdChanges(false);
		gallery.data.setFileName(to);
		if (blockto.contains("@"))
			gallery.data.selectBlock(blockto.substring(0, blockto.lastIndexOf("@")));
		this.saved = true;
//		reloadFile(to, false);
		reloadCombos();
		setGalleryTitle();
	}

	private void save() throws Exception
	{
		if (!saved)
			saveAs();
		else
		{
			if (dlgSave.saveActiveMetadataOnly())
				saveMd();
			else
				saveAll();
		}
	}// function save

	private void saveAs() throws Exception
	{
		if (dlgSave == null)
			dlgSave = new SaveJDialog(this, data.getMdSaveFileName(), false);
		else
			dlgSave.setMdFilename(data.getMdFilename());
		boolean save = dlgSave.showDialog(); // displays dialog and waits until
												// save or cancel clicked
		if (save)
		{
			if (dlgSave.saveActiveMetadataOnly())
				saveMd();
			else
				saveAll();

			if (dlgSave.doSaveImages())
				data.md.writeImages(dlgSave.getOutput(), dlgSave.isOutputIndependent(), dlgSave.getImageLabel());
		}

	}

	public void openMicrographs()
	{
		if (extractframe == null || !extractframe.isVisible())
			extractframe = ExtractParticlePicker.open(getBlock(), data.getFileName(), gallery.getImageWidth(), this);
		refreshExtractFrame();
	}

	private void refreshExtractFrame()
	{

		if (extractframe == null || !extractframe.isVisible())
			return;
		for (int i = 0; i < data.selection.length; i++)
			if (data.selection[i])
				extractframe.refreshActive(data.ids[i], data.isEnabled(i));

	}

	public void refreshActive(long id, boolean enabled)
	{
		for (int i = 0; i < data.ids.length; i++)
			if (id == data.ids[i])
			{
				goToImage(i);
				gallery.setSelectionEnabled(enabled);
			}

	}

	public String getBlock()
	{
		return (String) jcbBlocks.getSelectedItem();
	}

	public void reloadFile(String file, boolean changed) throws Exception
	{
		createModel();
		reloadMd(changed);
		reloadCombos();
	}
	
	private void reloadCombos(){
		createCombos();
		updateVisibleCombos();
		jcbBlocks.setSelectedItem(gallery.data.selectedBlock);
		
	}

	protected void initResliceButtonMenu()
	{
		// Create the popup menu.
		String reslice;
		final JPopupMenu popup = new JPopupMenu();
		JRadioButtonMenuItem mi;
		reslicegroup = new ButtonGroup();
		for (int i = 0; i < reslices.length; i++)
		{
			reslice = reslices[i];
			mi = new JRadioButtonMenuItem(reslice);
			reslicegroup.add(mi);
			if (i == 0)
				reslicegroup.setSelected(mi.getModel(), true);
			mi.setActionCommand(String.valueOf(ImageGeneric.VIEWS[i]));
			popup.add(mi);
			mi.addActionListener(new ResliceActionListener());

		}
		reslicebt = new JButton(XmippResource.getIcon("topview.png"));
		reslicebt.setToolTipText("Reslice");
		reslicebt.setContentAreaFilled(false);
		reslicebt.setFocusPainted(false);
		reslicebt.setOpaque(false);

		reslicebt.addMouseListener(new MouseAdapter()
		{
			public void mousePressed(MouseEvent e)
			{
				popup.show(e.getComponent(), e.getX(), e.getY());
			}
		});

	}
	
	
		
	

	class ResliceActionListener implements ActionListener
	{

		@Override
		public void actionPerformed(ActionEvent e)
		{
			JRadioButtonMenuItem mi = (JRadioButtonMenuItem) e.getSource();
			int view = Integer.parseInt(mi.getActionCommand());
			setResliceView(view);
			reslicegroup.setSelected(mi.getModel(), true);
		}

	}

	public Map<String, String> getKeyAssist()
	{
		Map<String, String> map = Collections.synchronizedMap(new LinkedHashMap<String, String>());
		map.put("Shift + Scroll Up/Ctrl + P", "Zoom in if images displayed");
		map.put("Shift + Scroll Down/Ctrl + O", "Zoom out if images displayed");
		map.put("Left click", "Selects a cell in gallery mode and a row in table mode");
		map.put("Right click", "Selects a row in table mode and displays row menu");
		map.put("Supr", "Delete selected cell in gallery mode and row in table mode");
		map.put("Up", "Select  previous row in table mode and cell in previous row in gallery mode");
		map.put("Down", "Select next row in table mode and cell in next row in gallery mode");
		map.put("Left", "Select  previous cell in gallery mode");
		map.put("Right", "Select next cell in gallery mode");
		map.put("Ctrl + O", "Open file");
		map.put("Ctrl + J", "Open with ImageJ");
		map.put("Ctrl + S", "Save");
		map.put("F5", "Refresh");
		map.put("Ctrl + Q", "Exit");
		map.put("Ctrl + N", "Display global normalization");
		map.put("Ctrl + L", "Display labels");
		map.put("Ctrl + R", "Render images");
		map.put("Ctrl + G", "Apply geometry");
		map.put("Ctrl + W", "Wrap");

		return map;
	}
        
       public boolean isImageSelected()
       {
           if(!data.allowGallery())
               return false;
           //what if there are no images on metadata??
           return gallery.getSelectionCount() > 0;
       }

}// class JFrameGallery
