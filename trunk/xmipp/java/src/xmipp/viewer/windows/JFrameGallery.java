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
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.File;
import java.util.ArrayList;

import javax.swing.AbstractAction;
import javax.swing.AbstractButton;
import javax.swing.ActionMap;
import javax.swing.ComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import xmipp.utils.XmippFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
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

import xmipp.ij.commons.ImagePlusLoader;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.Param;
import xmipp.utils.XmippWindowUtil;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippLabel;
import xmipp.utils.XmippMenuBarCreator;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.utils.XmippResource;
import xmipp.viewer.RowHeaderRenderer;
import xmipp.viewer.ctf.TasksEngine;
import xmipp.viewer.ctf.iCTFGUI;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.models.GalleryData;
import xmipp.viewer.models.GalleryRowHeaderModel;
import xmipp.viewer.models.ImageGallery;
import xmipp.viewer.models.MetadataGallery;
import xmipp.viewer.models.MicrographsTable;
import xmipp.viewer.windows.ClassesJDialog;

public class JFrameGallery extends JFrame implements iCTFGUI, WindowListener {
	private static final long serialVersionUID = -8957336972082018823L;

	private final static int DELAY_TO_UPDATE = 500;
	private static int update_counter = 0;
	public ImageGallery gallery;
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
	private JLabel jlGoto;
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

	protected static final float MAX_HEIGHT_RATE = 2.0f / 3.0f;
	// this rate is width/height
	protected static final float DIM_RATE = 4.0f / 3.0f;
	protected static final int MIN_WIDTH = 600;
	protected static int MIN_HEIGHT;
	protected static int MAX_HEIGHT;
	protected static int MAX_WIDTH;

	/** Some static initialization for fancy default dimensions */
	static {
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		float aux = (float) screenSize.height * MAX_HEIGHT_RATE;
		MAX_HEIGHT = Math.round(aux);
		aux = (float) MIN_WIDTH / DIM_RATE;
		MIN_HEIGHT = Math.round(aux);
		aux = (float) MAX_HEIGHT * DIM_RATE;
		MAX_WIDTH = Math.round(aux);
	}
	/** Store data about visualization */
	GalleryData data;

	/** Initialization function after GalleryData structure is created */
	private void init(GalleryData data) {
		try {
			this.data = data;
			createModel();
			createGUI();
		} catch (Exception e) {
			DEBUG.printException(e);
		}
	}

	/** Constructors */
	public JFrameGallery(String filename, Param parameters) {
		super();
		init(new GalleryData(filename, parameters, null));
	}

	public JFrameGallery(String filename, MetaData md, Param parameters) {
		super();
		init(new GalleryData(filename, parameters, md));
	}

	public JFrameGallery(String filenames[], Param parameters) {
		this(filenames, null, parameters);
	}

	public JFrameGallery(String filenames[], boolean enabled[], Param parameters) {
		super();
		// createGUI(new MDTableModel(filenames, enabled), parameters);
	}

	/**
	 * Open another metadata separataly *
	 */
	public void openMetadata(MetaData md) {
		new JFrameGallery(null, md, new Param());
	}

	/**
	 * Function to create the gallery type depending on the filename
	 * 
	 * @throws Exception
	 */
	private void createModel() throws Exception {
		gallery = data.createModel();
	}

	public GalleryData getData() {
		return data;
	}

	/**
	 * Function to create general GUI base on a TableModel. It will use helper
	 * functions to create different components of the GUI
	 */
	private void createGUI() {
		// Create file chooser and set current dir
		setIconImage(XmippResource.getIcon("xmipp_logo.png").getImage());
		fc = new XmippFileChooser();
		ctfTasks = new TasksEngine(JFrameGallery.this);

		isUpdating = true; // avoid handling some changes events

		setTitle(gallery.getTitle());
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setMinimumSize(new Dimension(MIN_WIDTH, MIN_HEIGHT));
		// addWindowListener(this);

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
		updateCombos();

		jspContent = new GalleryScroll();
		// Create table
		createTable();
		c.fill = GridBagConstraints.BOTH;
		c.gridx = 0;
		c.gridy = 2;
		c.weightx = 1.0;
		c.weighty = 1.0;
		container.add(jspContent, c);

		// Create the menu for table
		menu = new GalleryMenu();
		setJMenuBar(menu.getMenuBar());
		jpopUpMenuTable = new GalleryPopupMenu();
		menu.update();

		// pack();
		isUpdating = false;

		// Zoom in with Ctrl + P
		InputMap imap = container
				.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
		ActionMap amap = container.getActionMap();
		imap.put(KeyStroke.getKeyStroke("ctrl released P"), "zoomIn");
		amap.put("zoomIn", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				zoomChange(true);
			}

		});
		// Zoom in with Ctrl + O
		imap.put(KeyStroke.getKeyStroke("ctrl released M"), "zoomOut");
		amap.put("zoomOut", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				zoomChange(false);
			}

		});

		// Change view with Ctrl + Tab
		imap.put(KeyStroke.getKeyStroke("ctrl released I"), "changeView");
		amap.put("changeView", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				data.changeMode();
				reloadTableData();
			}
		});

		pack();
		XmippWindowUtil.centerWindows(this);
		setVisible(true);
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				addComponentListener(new java.awt.event.ComponentAdapter() {
					public void componentResized(
							java.awt.event.ComponentEvent evt) {
						formComponentResized(evt);
					}
				});
			}
		});
	}

	private void setInitialValues() {
		boolean adjust = false;
		if (data.parameters.columns > 0)
			gallery.setColumns(data.parameters.columns);
		else if (data.parameters.rows > 0)
			gallery.setRows(data.parameters.rows);
		else if (!data.isRotSpectraMode())
			adjust = true;
		if (data.isMicrographsMode()) {
			setExtendedState(JFrame.MAXIMIZED_BOTH);
			width = getSize().width;
		} else {
			int desiredCols = adjust ? (int) Math.ceil(Math.sqrt(gallery
					.getSize())) : gallery.getColumnCount();
			int desiredWidth = desiredCols * gallery.cellDim.width + 50;
			width = Math.min(Math.max(desiredWidth, MIN_WIDTH), MAX_WIDTH);
		}
		setAutoAdjustColumns(adjust);
	}

	/** Some tweaks over traditional JTable */
	public class GalleryScroll extends JScrollPane {
	}// class GalleryTable

	private void createTable() {
		// Create row header for enumerate rows
		try {
			rowHeaderModel = data.md.isColumnFormat() ? new GalleryRowHeaderModel(
					gallery.getRowCount(), 1) : new GalleryRowHeaderModel(data);
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		rowHeader = new JList();
		rowHeader.setModel(rowHeaderModel);
		LookAndFeel.installColorsAndFont(rowHeader, "TableHeader.background",
				"TableHeader.foreground", "TableHeader.font");
		rowHeader.setCellRenderer(new RowHeaderRenderer());
		jspContent.setRowHeaderView(rowHeader);

		table = new JTable();
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
		gallery.addTableModelListener(new TableModelListener() {

			@Override
			public void tableChanged(TableModelEvent e) {
				updateTable();
			}
		});

		table.addMouseListener(new java.awt.event.MouseAdapter() {
			public void mouseClicked(java.awt.event.MouseEvent evt) {
				tableMouseClicked(evt);
			}
		});

		// Zoom with Ctrl + MouseWeel
		table.addMouseWheelListener(new MouseWheelListener() {
			@Override
			public void mouseWheelMoved(MouseWheelEvent evt) {
				if (evt.isControlDown())
					zoomChange(evt.getWheelRotation() < 0);
				else
					table.getParent().dispatchEvent(evt);
			}
		});

		updateViewState();
		if (!adjustColumns())
			updateTable(); // update table if columns have not changed
	}// function createTable

	private void zoomChange(boolean increase) {
		int deltha = increase ? 10 : -10;
		jsZoom.setValue((Integer) jsZoom.getValue() + deltha);
	}

	private void updateTable() {
		// if (table.isShowing()) {
		boolean updatingState = isUpdating;
		isUpdating = true;
		update_counter++;
		DEBUG.printMessage(" *** Updating table: " + update_counter); // );
		// DEBUG.printStackTrace();

		// FIXME:gallery.updateSort();

		if (gallery.getSize() > 0) {
			// DEBUG.printMessage(String.format("updateTable: Table Model:\nsize: %d, cols: %d rows: %d",
			// gallery.getSize(), gallery.getColumnCount(),
			// gallery.getRowCount()));
			Dimension dimension = gallery.getCellSize();
			// renderer.setPreferredSize(dimension);

			// Adjusts rows size.
			table.setRowHeight(dimension.height);
			rowHeader.setFixedCellHeight(dimension.height);
			rowHeaderModel.setSize(gallery.getRowCount());
			// Adjusts columns width
			gallery.getColumnModel().adjustColumnsWidth(table);
			// columnModel.setWidth(dimension.width);

			// If auto adjust columns is enabled, refresh!
			jsRows.setValue(gallery.getRowCount());
			jsColumns.setValue(gallery.getColumnCount());

			rowHeader.revalidate();
			rowHeader.repaint();
			// table.revalidate();
			// repaint();
		}

		jsZoom.setValue(data.zoom);
		isUpdating = updatingState;
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				gallery.updateTableSelection(table);
			}
		});

		// }
	}// function updateTable

	/** Adjust the columns depending on the current windows width and cell width */
	private boolean adjustColumns() {
		if (autoAdjustColumns)
			return gallery.adjustColumn(width - 50);
		return false;
		// DEBUG.printMessage(String.format(
		// "==>> JFrameGallery.autoAdjust: width: %d", width));
		// int rw = rowHeader.getWidth();
		// DEBUG.printStackTrace();
		// FIXME
		// gallery.autoAdjustColumns(
		// // jsPanel.getVisibleRect().width - rowHeader.getWidth(),
		// // jsPanel.getViewportBorderBounds().width -
		// // rowHeader.getWidth(),
		// // jspContent.getViewport().getWidth() - rowHeader.getWidth()
		// w - rw, table.getIntercellSpacing().width);
		// updateColumnsRowsValues();
	}

	public void setAutoAdjustColumns(boolean autoAdjustColumns) {
		this.autoAdjustColumns = autoAdjustColumns;
		jcbAutoAdjustColumns.setSelected(autoAdjustColumns);
		jsColumns.setEnabled(!autoAdjustColumns);
		jsRows.setEnabled(!autoAdjustColumns);
	}

	private void goToImage(int index) {
		gallery.gotoItem(index);

		int coords[] = gallery.getCoords(index);
		DEBUG.printMessage(String.format(
				"gotoImage, index: %d, row: %d, col:%d", index, coords[0],
				coords[1]));

		// Gets current selected cell bounds.
		Rectangle rect = table.getCellRect(coords[0], coords[1], true);

		// Ensures item is visible
		Point pos = jspContent.getViewport().getViewPosition();
		rect.translate(-pos.x, -pos.y);
		jspContent.getViewport().scrollRectToVisible(rect);

		repaint();
	}

	public class Worker implements Runnable {
		public static final int STATS = 0;
		public static final int PCA = 1;
		public static final int FSC = 2;
		public String message;
		/** Constructor selecting operation */
		private int op; // store operation

		public Worker(int operation) {
			op = operation;
		}

		public void run() {
			try {
				switch (op) {
				case STATS:
					computeStatsImages();
					break;
				case PCA:
					pca();
					break;
				case FSC:
					fsc();
					break;
				}
			} catch (Exception e) {
				showException(e);
			}
			ImagesWindowFactory.releaseGUI(JFrameGallery.this.getRootPane());
		}

		public String getMessage() {
			switch (op) {
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
	public void runInBackground(int operation) {
		Worker w = new Worker(operation);
		ImagesWindowFactory.blockGUI(getRootPane(), w.getMessage());
		Thread thr = new Thread(w);
		thr.start();
	}

	private void computeStatsImages() throws Exception {
		ImageGeneric imgAvg = new ImageGeneric();
		ImageGeneric imgStd = new ImageGeneric();
		data.md.getStatsImages(imgAvg, imgStd, data.useGeo,
				data.getRenderLabel());
		ImagePlus impAvg = XmippImageConverter.convertToImagePlus(imgAvg);
		ImagePlus impStd = XmippImageConverter.convertToImagePlus(imgStd);
		imgAvg.destroy();
		imgStd.destroy();
		XmippImageWindow winAvg = new XmippImageWindow(new ImagePlusLoader(impAvg), "AVG: "
				+ data.filename);
		XmippWindowUtil.setLocation(0.2f, 0.5f, winAvg, this);
		winAvg.setVisible(true);
		XmippImageWindow winStd = new XmippImageWindow(new ImagePlusLoader(impStd), "STD: "
				+ data.filename);
		XmippWindowUtil.setLocation(0.8f, 0.5f, winStd, this);
		winStd.setVisible(true);
	}

	private void save() {
		if (dlgSave == null)
			dlgSave = new SaveJDialog(this);
		if (!saved)
			saveAs();
		else {
			try {
				String path = dlgSave.getMdFilename();
				if (dlgSave.isAppendMode())
					data.md.writeBlock(path);
				else
					data.md.write(path);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	private void saveAs() {
		if (dlgSave == null)
			dlgSave = new SaveJDialog(this);
		dlgSave.setMdFilename(data.getMdFilename());
		if (dlgSave.showDialog()) {
			try {
				String path = dlgSave.getMdFilename();
				DEBUG.printMessage(path);
				if (dlgSave.isAppendMode())
					data.md.writeBlock(path);
				else
					data.md.write(path);
				saved = true;
				if (dlgSave.doSaveImages()) {
					data.md.writeImages(dlgSave.getOutput(),
							dlgSave.isOutputIndependent(),
							dlgSave.getImageLabel());
				}
			} catch (Exception e) {
				XmippDialog.showError(this, e.getMessage());
			}
		}
	}

	private boolean openClassesDialog() {
		if (dlgClasses == null) {
			dlgClasses = new ClassesJDialog(JFrameGallery.this);
		}
		boolean result = dlgClasses.showDialog();
		dlgClasses.resetClasses();
		return result;
	}

	static String forceExtension(String filename, String ext) {
		int dot = filename.lastIndexOf(".");
		return filename.substring(0, dot) + ext;
	}

	public void pca() throws Exception {
		ImageGeneric image = new ImageGeneric();
		data.md.getPCAbasis(image, data.getRenderLabel());
		ImagePlus imp = XmippImageConverter.convertToImagePlus(image);
		// new XmippImageWindow(imp, "PCA: " + data.filename);
		imp.setTitle("PCA: " + data.filename);
		ImagesWindowFactory.captureFrame(imp);

	}

	public void fsc() throws Exception {
		JFrameFSC frame = new JFrameFSC(data);
		XmippWindowUtil.centerWindows(frame, this);
		frame.setVisible(true);
	}

	/***
	 * Helper function to create toolbar toggle buttons
	 */
	protected void setupButton(AbstractButton btn, String icon, String text,
			ActionListener listener) {
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

	private void updateViewState() {
		ImageIcon icon;
		String text;
		if (!data.isGalleryMode()) {
			icon = XmippResource.VIEW_GALLERY_ICON;
			text = XmippLabel.LABEL_VIEW_GALLERY;
		} else {
			icon = XmippResource.VIEW_MD_ICON;
			text = XmippLabel.LABEL_VIEW_MD;
		}
		btnChangeView.setIcon(icon);
		btnChangeView.setToolTipText(text);
		if (data.isTableMode()) { // if we are in table mode only allow change
									// if exist render label
			boolean hasRender = data.allowGallery();
			btnChangeView.setEnabled(hasRender);
			jsZoom.setEnabled(hasRender);
			jlZoom.setEnabled(hasRender);
		}
	}

	/** Reload table data */
	public void reloadTableData() {
		try {
			table.removeAll();
			createModel();
			// gallery.setShowLabels(menu.getShowLabel());
			createTable();

			menu.update();
			updateCombos();
			if (dlgSave != null)
				dlgSave.setInitialValues();
			saved = false;
			setTitle(gallery.getTitle());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/** Reload metadata info, rebuild the table */
	private void reloadMd() throws Exception {
		data.loadMd();
		reloadTableData();
	}// function reloadMd

	/**
	 * Fill some label mode can be: "constant", "linear", "uniform", "gaussian"
	 * values is a list of string
	 * */
	public void fillLabel(int label, String mode, String... values)
			throws Exception {
		DEBUG.printFormat("fillLabels, mode:%s\n", mode);
		if (mode.equalsIgnoreCase(MetaData.FILL_CONSTANT))
			data.md.fillConstant(label, values[0]);
		else {
			Double v1 = Double.parseDouble(values[0]);
			Double v2 = Double.parseDouble(values[1]);
			if (mode.equalsIgnoreCase(MetaData.FILL_LINEAR))
				data.md.fillLinear(label, v1, v2);
			else if (mode.equalsIgnoreCase(MetaData.FILL_RAND_UNIFORM))
				data.md.fillRandom(label, "uniform", v1, v2);
			else if (mode.equalsIgnoreCase(MetaData.FILL_RAND_GAUSSIAN))
				data.md.fillRandom(label, "gaussian", v1, v2);
		}
		reloadMd();
	}

	/**
	 * Delete selected or disabled items if 'selected' is true, selection is
	 * removed if false, the disabled items
	 * */
	public void removeObjects(boolean selected) throws Exception {
		String type = selected ? "selected" : "disabled";
		if (XmippDialog.showWarning(this,
				String.format("Are you sure to delete %s items?", type))) {
			if (selected)
				data.removeSelection();
			else
				data.md.removeDisabled();
			reloadMd();
		}
	}

	/** Drop some label from the metadata */
	public void removeLabel(int label) throws Exception {
		data.md.removeLabel(label);
		reloadMd();
	}

	/***
	 * Function to create the main toolbar
	 */
	protected void createToolbar() {
		// Create Main TOOLBAR
		toolBar = new JToolBar();
		toolBar.setRollover(true);
		toolBar.setLayout(new FlowLayout(FlowLayout.LEFT));

		btnChangeView = new JButton();
		// updateViewState();
		btnChangeView.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				data.changeMode();
				reloadTableData();
			}
		});

		toolBar.add(btnChangeView);
		toolBar.addSeparator();

		jlZoom = new javax.swing.JLabel();
		jsZoom = new javax.swing.JSpinner();
		jlZoom.setIcon(XmippResource.getIcon(XmippResource.ZOOM));
		jlZoom.setToolTipText(XmippLabel.LABEL_ZOOM);
		toolBar.add(jlZoom);

		jsZoom.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(1),
				Integer.valueOf(1), null, Integer.valueOf(1)));
		jsZoom.addChangeListener(new javax.swing.event.ChangeListener() {
			public void stateChanged(javax.swing.event.ChangeEvent evt) {
				Integer zoom = (Integer) jsZoom.getValue();
				gallery.setZoom(zoom);
				// gallery.updateTableSelection(table);
			}
		});

		toolBar.add(jsZoom);
		toolBar.addSeparator();

		jlGoto = new javax.swing.JLabel();
		jsGoToImage = new javax.swing.JSpinner();
		jlGoto.setIcon(XmippResource.getIcon(XmippResource.GOTO));
		jlGoto.setToolTipText(XmippLabel.LABEL_GOTO_ITEM);
		toolBar.add(jlGoto);

		jsGoToImage.setValue(1);
		jsGoToImage.addChangeListener(new javax.swing.event.ChangeListener() {
			public void stateChanged(javax.swing.event.ChangeEvent evt) {
				jsGoToImageStateChanged(evt);
			}
		});
		toolBar.add(jsGoToImage);

		toolBar.addSeparator();

		jcbAutoAdjustColumns = new JToggleButton();
		setupButton(jcbAutoAdjustColumns, XmippResource.ADJUST_COLS,
				XmippLabel.MSG_ADJUST_COLS,
				new java.awt.event.ActionListener() {
					public void actionPerformed(java.awt.event.ActionEvent evt) {
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

		jsColumns.addChangeListener(new javax.swing.event.ChangeListener() {
			public void stateChanged(javax.swing.event.ChangeEvent evt) {
				jsColumnsStateChanged(evt);
			}
		});
		toolBar.add(jsColumns);

		jlRows.setText(XmippLabel.LABEL_ROWS);
		toolBar.add(jlRows);

		jsRows.addChangeListener(new javax.swing.event.ChangeListener() {
			public void stateChanged(javax.swing.event.ChangeEvent evt) {
				jsRowsStateChanged(evt);
			}
		});
		toolBar.add(jsRows);

		// Some settings of the spinners
		jsRows.setModel(new SpinnerNumberModel(1, 1, gallery.getSize(), 1));
		jsColumns.setModel(new SpinnerNumberModel(1, 1, gallery.getSize(), 1));
		jsGoToImage
				.setModel(new SpinnerNumberModel(1, 1, gallery.getSize(), 1));

		int TEXTWIDTH = 4;
		((JSpinner.NumberEditor) jsZoom.getEditor()).getTextField().setColumns(
				TEXTWIDTH);
		((JSpinner.NumberEditor) jsGoToImage.getEditor()).getTextField()
				.setColumns(TEXTWIDTH);
		((JSpinner.NumberEditor) jsRows.getEditor()).getTextField().setColumns(
				TEXTWIDTH);
		((JSpinner.NumberEditor) jsColumns.getEditor()).getTextField()
				.setColumns(TEXTWIDTH);

	}// function createToolbar

	/** Create combos for selection of block and volume if its the case */
	protected void createCombos() {
		cbPanel = new JPanel();
		cbPanel.setLayout(new FlowLayout(FlowLayout.LEFT));

		// Add blocks selector combo
		jlBlocks = new JLabel(XmippLabel.LABEL_BLOCK);
		cbPanel.add(jlBlocks);
		jcbBlocks = new JComboBox();
		jcbBlocks.setModel(new ComboBoxModel() {

			@Override
			public int getSize() {
				return data.mdBlocks.length;
			}

			@Override
			public Object getElementAt(int index) {
				return data.mdBlocks[index];
			}

			@Override
			public void setSelectedItem(Object item) {
				data.selectBlock((String) item);
				reloadTableData();
			}

			@Override
			public Object getSelectedItem() {
				return data.selectedBlock;
			}

			@Override
			public void removeListDataListener(ListDataListener arg0) {
			}

			@Override
			public void addListDataListener(ListDataListener arg0) {
				// TODO Auto-generated method stub
			}
		});
		cbPanel.add(jcbBlocks);
		// Add volumes selector combo
		jlVolumes = new JLabel(XmippLabel.LABEL_VOLUME);
		cbPanel.add(jlVolumes);
		jcbVolumes = new JComboBox();
		jcbVolumes.setModel(new ComboBoxModel() {

			@Override
			public int getSize() {
				return data.getNumberOfVols();
			}

			@Override
			public Object getElementAt(int index) {
				return data.getVolumeAt(index);
			}

			@Override
			public void setSelectedItem(Object anItem) {
				data.selectVolume((String) anItem);
				reloadTableData();
			}

			@Override
			public Object getSelectedItem() {
				return data.selectedVolFn;
			}

			@Override
			public void addListDataListener(ListDataListener arg0) {
				// TODO Auto-generated method stub

			}

			@Override
			public void removeListDataListener(ListDataListener arg0) {
				// TODO Auto-generated method stub

			}
		});
		cbPanel.add(jcbVolumes);
	}

	protected void updateCombos() {
		boolean showBlocks = data.getNumberOfBlocks() > 1;
		boolean showVols = data.getNumberOfVols() > 1 && data.isVolumeMode();
		jcbBlocks.setVisible(showBlocks);
		jcbVolumes.setVisible(showVols);
		jlBlocks.setVisible(showBlocks);
		jlVolumes.setVisible(showVols);
		cbPanel.setVisible(showBlocks || showVols);

	}

	private void jsRowsStateChanged(javax.swing.event.ChangeEvent evt) {// GEN-FIRST:event_jsRowsStateChanged
		if (!isUpdating) {
			gallery.setRows((Integer) jsRows.getValue());
		}
	}

	private void jsColumnsStateChanged(javax.swing.event.ChangeEvent evt) {// GEN-FIRST:event_jsColumnsStateChanged
		if (!isUpdating) {
			gallery.setColumns((Integer) jsColumns.getValue());
		}
	}

	private void jsGoToImageStateChanged(javax.swing.event.ChangeEvent evt) {
		DEBUG.printMessage("jsGotoImage...");
		if (!isUpdating) {
			goToImage((Integer) jsGoToImage.getValue() - 1);
		}
	}

	private void formComponentResized(java.awt.event.ComponentEvent evt) {
		width = getSize().width;
		if (!isUpdating && autoAdjustColumns)
			adjustColumns();
	}

	private void tableMouseClicked(MouseEvent evt) {
		final Point p = evt.getPoint();
		int view_row = table.rowAtPoint(p);
		int view_col = table.columnAtPoint(p);

		if (evt.getButton() == MouseEvent.BUTTON1) { // Left click.
			if (evt.getClickCount() > 1) {
				try {
					// String imageFn = gallery.getImageFilenameAt(view_row,
					// view_col);
					// if (imageFn != null)
					// new XmippImageWindow(imageFn);
					gallery.handleDoubleClick(view_row, view_col);
				} catch (Exception e) {
					XmippDialog.showError(this, e.getMessage());
				}
			} else {
				// Ctrl adds items to selection, otherwise previous ones are
				// removed.
				boolean move = true;
				if (!evt.isControlDown() && !evt.isShiftDown()) {
					if (gallery.getSelectionCount() <= 1
							|| XmippDialog
									.showWarning(this,
											"You will lose previous selection.\nDo you want to proceed?")) {
						gallery.clearSelection();
						gallery.touchItem(view_row, view_col);

					} else
						move = false;
				}

				else {

					if (evt.isShiftDown()) {
						gallery.selectRange(previousSelectedRow,
								previousSelectedCol, view_row, view_col, true);
					} else if (evt.isControlDown()) {
						gallery.touchItem(view_row, view_col);
					}
				}
				if (move) {
					isUpdating = true;
					jsGoToImage
							.setValue(gallery.getIndex(view_row, view_col) + 1);
					isUpdating = false;
				}
				if (!evt.isShiftDown()) {
					previousSelectedRow = view_row;
					previousSelectedCol = view_col;
				}
			}

		} else if (evt.getButton() == MouseEvent.BUTTON3) { // Right click.

			final MouseEvent me = evt;
			if (gallery.handleRightClick(view_row, view_col, jpopUpMenuTable)) {
				if (gallery.getSelectionCount() < 2) {
					gallery.clearSelection();
					gallery.touchItem(view_row, view_col);
				}
				SwingUtilities.invokeLater(new Runnable() {
					public void run() {
						jpopUpMenuTable.show(me.getComponent(), p);
					}
				});
			}
		}
		table.invalidate();
		table.repaint();
	}// function tableMouseClicked

	private void autoAdjustColumns(boolean value) {
		setAutoAdjustColumns(value);
		adjustColumns();
	}
	
	private void setResliceView(int view){
		data.resliceView = view;
		reloadTableData();
	}

	class GalleryMenu extends XmippMenuBarCreator {

		@Override
		protected void createItems() throws Exception {
			// File
			addItem(FILE, "File");
			addItem(FILE_OPEN, "Open ...", null, "control released O");
			addItem(FILE_OPENWITH_IJ, "Open with ImageJ", "ij.gif",
					"control released J");
			addItem(FILE_OPENWITH_CHIMERA, "Open with Chimera", "chimera.gif",
					"control released H");
			addSeparator(FILE);
			addItem(FILE_SAVE, "Save", "save.gif", "control released S");
			addItem(FILE_SAVEAS, "Save as", "save_as.gif");
			addItem(FILE_REFRESH, "Refresh", "refresh.gif", "released F5");
			addSeparator(FILE);
			addItem(FILE_EXIT, "Exit", null, "control released Q");
			// Display
			addItem(DISPLAY, "Display");
			addItem(DISPLAY_NORMALIZE, "Global normalization", null,
					"control released N");
			addItem(DISPLAY_SHOWLABELS, "Show labels", null,
					"control released L");
			addSeparator(DISPLAY);
			addItem(DISPLAY_RENDERIMAGES, "Render images", null,
					"control released R");
			addItem(DISPLAY_APPLYGEO, "Apply geometry", null,
					"control released G");
			addItem(DISPLAY_WRAP, "Wrap", null, "control released W");
			addItem(DISPLAY_COLUMNS, "Columns ...", "columns.gif");
			addItem(DISPLAY_RESLICE, "Reslice");
			String text[] = {"Z Negative (Front)", "Y Negative (Top)", "X Negative (Left)", "Y Positive (Bottom)", "X Positive (Right)" };
			for (int i = 0; i < ImageGeneric.VIEWS.length; ++i)
				addItem(DISPLAY_RESLICE_VIEWS[i], text[i]);
			// Metadata operations
			addItem(METADATA, "Metadata");
			addItem(STATS, "Statistics");
			addItem(STATS_AVGSTD, "Avg & Std images");
			addItem(STATS_PCA, "PCA");
			addItem(STATS_FSC, "FSC");
			addItem(MD_CLASSES, "Superclasses");
			addItem(MD_EDIT_COLS, "Edit labels", "edit.gif");
			addItem(MD_ADD_OBJECT, "Add new object", "new_object.gif");
			addItem(MD_REMOVE_DISABLED, "Remove disabled", "delete.gif");
			addItem(MD_REMOVE_SELECTION, "Remove selection");
			// Help
			addItem(HELP, "Help");
			addItem(HELP_ONLINE, "Online help", "online_help.gif");
		}// function createItems

		public void update() {
			boolean galMode = data.isGalleryMode();
			boolean volMode = data.isVolumeMode();
			setItemEnabled(FILE_OPENWITH_CHIMERA, volMode);
			setItemEnabled(FILE_SAVE, !volMode);
			setItemEnabled(FILE_SAVEAS, !volMode);
			setItemSelected(DISPLAY_NORMALIZE, gallery.getNormalized());
			setItemEnabled(DISPLAY_APPLYGEO, data.containsGeometryInfo());
			setItemEnabled(DISPLAY_WRAP, data.containsGeometryInfo()
					&& data.useGeo);
			setItemSelected(DISPLAY_WRAP, data.wrap);
			setItemSelected(DISPLAY_APPLYGEO, data.useGeo);
			setItemSelected(DISPLAY_APPLYGEO, data.useGeo && data.wrap);
			setItemEnabled(DISPLAY_RENDERIMAGES,
					!galMode && data.hasRenderLabel());
			setItemSelected(DISPLAY_RENDERIMAGES, data.globalRender);
			for (int i = 0; i < ImageGeneric.VIEWS.length; ++i)
				setItemSelected(DISPLAY_RESLICE_VIEWS[i], 
						(data.resliceView == ImageGeneric.VIEWS[i]));
			setItemEnabled(DISPLAY_COLUMNS, !galMode);
			setItemEnabled(DISPLAY_RESLICE, volMode);
			setItemEnabled(MD_CLASSES, data.is2DClassificationMd());
			setItemEnabled(MD_REMOVE_SELECTION, gallery.getSelectionCount() > 0);
			setItemEnabled(STATS, !volMode);
		}// function update

		@Override
		protected void handleActionPerformed(ActionEvent evt) {
			String cmd = evt.getActionCommand();
			try {
				if (cmd.equals(DISPLAY_NORMALIZE)) {
					gallery.setNormalized(getItemSelected(DISPLAY_NORMALIZE));
				} else if (cmd.equals(DISPLAY_APPLYGEO)
						|| cmd.equals(DISPLAY_WRAP)) {
					if (data.containsGeometryInfo()) {
						((MetadataGallery) gallery).setUseGeometry(
								getItemSelected(DISPLAY_APPLYGEO),
								getItemSelected(DISPLAY_WRAP));
						setItemEnabled(DISPLAY_WRAP,
								data.containsGeometryInfo() && data.useGeo);
					}
				} else if (cmd.equals(DISPLAY_SHOWLABELS)) {
					gallery.setShowLabels(getItemSelected(DISPLAY_SHOWLABELS));
				} else if (cmd.equals(DISPLAY_RENDERIMAGES)) {
					gallery.setRenderImages(getItemSelected(DISPLAY_RENDERIMAGES));
				} else if (cmd.equals(DISPLAY_COLUMNS)) {
					ColumnsJDialog dialog = new ColumnsJDialog(
							JFrameGallery.this);
					boolean result = dialog.showDialog();
					if (result) {
						ArrayList<ColumnInfo> columns = dialog
								.getColumnsResult();
						isUpdating = true;
						((MetadataGallery) gallery).updateColumnInfo(columns);
						gallery.fireTableDataChanged();
						setItemEnabled(DISPLAY_RENDERIMAGES, data.globalRender);
						// menu.enableRenderImages(data.globalRender);
						isUpdating = false;
					}
				} else if (cmd.equals(STATS_AVGSTD))
					runInBackground(Worker.STATS);
				else if (cmd.equals(STATS_PCA))
					runInBackground(Worker.PCA);
				else if (cmd.equals(STATS_FSC))
					runInBackground(Worker.FSC);
				else if (cmd.equals(FILE_OPEN)) {
					if (fc.showOpenDialog(JFrameGallery.this) != XmippFileChooser.CANCEL_OPTION) {
						if (fc.getSelectedFile().exists())
							ImagesWindowFactory.openFileAsDefault(fc
									.getSelectedPath());
						else
							XmippDialog.showError(JFrameGallery.this, String
									.format("File: '%s' doesn't exist.",
											fc.getSelectedPath()));
					}
				} else if (cmd.equals(FILE_SAVE)) {
					save();
				} else if (cmd.equals(FILE_SAVEAS)) {
					saveAs();
				} else if (cmd.equals(FILE_EXIT)) {
					System.exit(0);
				} else if (cmd.equals(FILE_OPENWITH_CHIMERA)) {
					try {
						String args = data.selectedVolFn;
						if (Filename.isSpiderVolume(args))
							args = "spider:" + args;
						// FIXME: Check chimera is installed
						Process p = new ProcessBuilder("chimera", args).start();
					} catch (Exception ex) {
						ex.printStackTrace();
					}
				} else if (cmd.equals(FILE_OPENWITH_IJ)) {
					try {
						ImagePlus imp = gallery.getImagePlus();
						ImagesWindowFactory.captureFrame(imp);
					} catch (Exception e1) {
						e1.printStackTrace();
					}
				} else if (cmd.equals(FILE_REFRESH)) {
					data.reloadMd();
					reloadTableData();
				} else if (cmd.contains(DISPLAY_RESLICE)) {
					for (int i = 0; i < ImageGeneric.VIEWS.length; ++i)
						if (cmd.equals(DISPLAY_RESLICE_VIEWS[i])){
							setResliceView(ImageGeneric.VIEWS[i]);
							break;
						}
				} else if (cmd.equals(MD_CLASSES)) {
					openClassesDialog();
				} else if (cmd.equals(MD_EDIT_COLS)) {
					EditLabelsJDialog dlg = new EditLabelsJDialog(
							JFrameGallery.this);
					dlg.showDialog();
				} else if (cmd.equals(MD_REMOVE_SELECTION)) {
					removeObjects(true);
				} else if (cmd.equals(MD_REMOVE_DISABLED)) {
					removeObjects(false);
				} else if (cmd.equals(MD_ADD_OBJECT)) {
					AddObjectJDialog dlg = new AddObjectJDialog(JFrameGallery.this);
					if (dlg.showDialog()){
						data.md.unionAll(dlg.md);
						reloadMd();
					}
				} else if (cmd.equals(HELP_ONLINE)) {
					XmippWindowUtil.openURI("http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/WebHome");
				} 
				
			} catch (Exception e) {
				showException(e);
			}
		}// function handleActionPerformed
	}// class GalleryMenu

	class GalleryPopupMenu extends XmippPopupMenuCreator {
		protected int row;
		protected int col;

		@Override
		protected void createItems() throws Exception {
			addItem(ENABLED, "Enable", "enable.gif");
			addItem(DISABLED, "Disable", "disable.gif");
			addItem(REFRESH, "Refresh", "refresh.gif");
			// addSeparator();
			addItem(OPEN, "Open");
			addItem(OPEN_ASTEXT, "Open as text");
			addItem(CTF_PROFILE, "Show CTF profile");
			addItem(CTF_RECALCULATE, "Recalculate CTF");
			addSeparator();
			addItem(SET_CLASS, "Set class");
			addItem(OPEN_IMAGES, "Open images");
			addItem(SELECT, "Select");
			addItem(SELECT_ALL, "All", null, "control released A");
			addItem(SELECT_TOHERE, "To here");
			addItem(SELECT_FROMHERE, "From here");
			initItems();
		}// function createItems

		public void show(Component cmpnt, Point location) {
			// This item visibility depends on current selection
			setItemVisible(OPEN_IMAGES,
					data.is2DClassificationMd()
							&& gallery.getSelectionCount() == 1);
			// Update menu items status depending on item.
			row = table.rowAtPoint(location);
			col = table.columnAtPoint(location);
			getPopupMenu().show(cmpnt, location.x, location.y);

		}// function show

		private void selectRange(int first, int last) {
			gallery.selectRange(first, last, true);
		}

		private void showCTF(boolean profile) {
			try {
				String ctfModel = data.md.getValueString(MDLabel.MDL_CTFMODEL,
						data.ids[row]);
				String displayFilename = data.md.getValueString(
						MDLabel.MDL_ASSOCIATED_IMAGE2, data.ids[row]);
				String psdFile = data.md.getValueString(MDLabel.MDL_PSD,
						data.ids[row]);

				ImageGeneric img = new ImageGeneric(displayFilename);
				ImagePlus imp = XmippImageConverter.readToImagePlus(img);

				if (profile)
					ImagesWindowFactory.openCTFWindow(imp, ctfModel, psdFile);
				else
					ImagesWindowFactory.openCTFImage(imp, ctfModel, psdFile,
							ctfTasks, data.md.getFilename(), row);

			} catch (Exception e) {
				XmippDialog.showError(JFrameGallery.this, e.getMessage());
			}
		}

		/** Set values to defaults */
		@Override
		public void initItems() {
			setItemVisible(OPEN, false);
			setItemVisible(OPEN_ASTEXT, false);
			setItemVisible(CTF_PROFILE, false);
			setItemVisible(CTF_RECALCULATE, false);
			setItemVisible(SET_CLASS, data.is2DClassificationMd());
		}

		@Override
		protected void handleActionPerformed(ActionEvent evt) {
			String cmd = evt.getActionCommand();
			if (cmd.equals(SELECT_ALL)) {
				selectRange(0, gallery.getSize() - 1);
			} else if (cmd.equals(SELECT_TOHERE)) {
				selectRange(0, gallery.getIndex(row, col));
			} else if (cmd.equals(SELECT_FROMHERE)) {
				selectRange(gallery.getIndex(row, col), gallery.getSize() - 1);
			} else if (cmd.equals(ENABLED)) {
				gallery.setSelectionEnabled(true);
				gallery.clearSelection();
			} else if (cmd.equals(DISABLED)) {
				gallery.setSelectionEnabled(false);
				gallery.clearSelection();
			} else if (cmd.equals(REFRESH)) {
				gallery.refreshAt(row, col);
			} else if (cmd.equals(OPEN)) {
				String file = data.getValueFromCol(row, col);
				ImagesWindowFactory.openFileAsDefault(file);
			} else if (cmd.equals(OPEN_ASTEXT)) {
				String file = gallery.getValueAt(row, col).toString();
				ImagesWindowFactory.openFileAsText(file, null);
			} else if (cmd.equals(CTF_PROFILE)) {
				showCTF(true);
			} else if (cmd.equals(CTF_RECALCULATE)) {
				showCTF(false);
			} else if (cmd.equals(SET_CLASS)) {
				if (openClassesDialog()) {
					int classNumber = dlgClasses.getSelectedClass();
					// DEBUG.printMessage(String.format("class: %d",
					// classNumber));
					gallery.setSelectionClass(classNumber);
				}
			} else if (cmd.equals(OPEN_IMAGES)) {
				int index = gallery.getIndex(row, col);
				MetaData md = data.getClassImages(index);
				if (md != null)
					openMetadata(md);
				else
					XmippDialog.showWarning(JFrameGallery.this, "This class has no images");
			}
			initItems();
		}

	}// class JPopUpMenuGallery

	public void showException(Exception e) {
		XmippDialog.showException(this, e);
	}

	@Override
	public void setRunning(boolean running) {
		// XmippDialog.showInfo(this, String.format("Calculating ctf"));
		// ImagesWindowFactory.blockGUI(getRootPane(), "Calculating CTF");
	}

	@Override
	public void setRowBusy(int row) {
		((MicrographsTable) gallery).setRowBusy(row);
	}

	@Override
	public void setRowIdle(int row) {
		((MicrographsTable) gallery).setRowIdle(row);
	}

	@Override
	public String getFilename() {
		return data.getMdFilename();
	}

	@Override
	public void done() {
		XmippDialog.showInfo(this, String.format("Calculating ctf: DONE"));
	}

	@Override
	public void windowActivated(WindowEvent arg0) {
		// TODO Auto-generated method stub

	}

	@Override
	public void windowClosed(WindowEvent arg0) {
		// TODO Auto-generated method stub

	}

	@Override
	public void windowClosing(WindowEvent arg0) {
		// DEBUG.printMessage("on window closing...");
		// if (dlgClasses != null) {
		// DEBUG.printMessage("disposing dialog");
		// dlgClasses.dispose();
		// }

	}

	@Override
	public void windowDeactivated(WindowEvent arg0) {
		// TODO Auto-generated method stub

	}

	@Override
	public void windowDeiconified(WindowEvent arg0) {
		// TODO Auto-generated method stub

	}

	@Override
	public void windowIconified(WindowEvent arg0) {
		// TODO Auto-generated method stub

	}

	@Override
	public void windowOpened(WindowEvent arg0) {
		// TODO Auto-generated method stub

	}
}// class JFrameGallery
