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
package xmipp.viewer;

import ij.IJ;
import ij.ImagePlus;

import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.ItemEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Timer;
import java.util.TimerTask;

import javax.swing.AbstractAction;
import javax.swing.AbstractButton;
import javax.swing.ActionMap;
import javax.swing.ButtonGroup;
import javax.swing.ComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
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

import xmipp.ij.XmippImageConverter;
import xmipp.ij.XmippImageWindow;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.Param;
import xmipp.utils.WindowUtil;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippLabel;
import xmipp.utils.XmippMenuCreator;
import xmipp.utils.XmippResource;
import xmipp.viewer.windows.ImagesWindowFactory;

public class JFrameGallery extends JFrame {
	private static final long serialVersionUID = -8957336972082018823L;

	private final static int DELAY_TO_UPDATE = 500;
	private static int update_counter = 0;
	private ImageGallery gallery;
	private GalleryRowHeaderModel rowHeaderModel;
	private int previousSelectedRow, previousSelectedCol;
	private JList rowHeader;
	private Timer updateTimer = new Timer(true); // Timer for zoom.
	private TableUpdater tableUpdaterTask; // Associated task for zoom timer.
	private boolean isUpdating;
	// this flag will be used to avoid firing properties change events
	// when the change is from our code and not external user interaction
	private boolean autoAdjustColumns = false;
	private GalleryPopupMenu jpopUpMenuTable;
	private GalleryMenu menu;
	private JFileChooser fc; 
	private SaveJDialog dlgSave = null;

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
	// private javax.swing.JToggleButton jtbNormalize;
	// private javax.swing.JToggleButton jtbUseGeometry;
	private javax.swing.JTable table;
	private javax.swing.JToolBar toolBar;

	/** Store data about visualization */
	GalleryData data;

	public enum RESLICE_MODE {

		TOP_Y, RIGHT_X
	}

	/** Initialization function after GalleryData structure is created */
	private void init(GalleryData data) {
		try {
			this.data = data;
			createModel();
			createGUI();
		} catch (Exception e) {
			DEBUG.printException(e);
			IJ.error(e.getMessage());
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
		//Create file chooser and set current dir
		fc = new JFileChooser();
		fc.setCurrentDirectory(new File(System.getProperty("user.dir")));
		
		isUpdating = true; // avoid handling some changes events

		setTitle(gallery.getTitle());
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setMinimumSize(new Dimension(600, 400));
		setPreferredSize(new Dimension(800, 600));
		ImagesWindowFactory.setConvenientSize(this);

		// Get main pane and set layout
		Container container = getContentPane();
		container.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();

		// Create toolbar buttons
		createToolbar();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridy = 0;
		container.add(toolBar, c);

		// Create combos for selection of blocks and/or volumes
		createCombos();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridy = 1;
		container.add(cbPanel, c);
		updateCombos();

		jspContent = new GalleryScroll();
		setInitialValues();
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

		pack();
		isUpdating = false;

		addComponentListener(new java.awt.event.ComponentAdapter() {
			public void componentResized(java.awt.event.ComponentEvent evt) {
				formComponentResized(evt);
			}
		});
		setVisible(true);
	}

	private void setInitialValues() {
		boolean adjust = false;
		if (data.parameters.columns > 0)
			gallery.setColumns(data.parameters.columns);
		else if (data.parameters.rows > 0)
			gallery.setRows(data.parameters.rows);
		else {
			adjustColumns();
			adjust = true;
		}
		setAutoAdjustColumns(adjust);
	}

	/** Some tweaks over traditional JTable */
	public class GalleryScroll extends JScrollPane {
	}// class GalleryTable

	private void createTable() {
		// Create row header for enumerate rows
		rowHeaderModel = new GalleryRowHeaderModel(gallery.getRowCount(), 1);
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

		//Zoom with Ctrl + MouseWeel
		table.addMouseWheelListener(new MouseWheelListener() {
			@Override
			public void mouseWheelMoved(MouseWheelEvent evt) {
				if (evt.isControlDown())
					zoomChange(evt.getWheelRotation() < 0);
				else
					table.getParent().dispatchEvent(evt);
			}
		});
		
		//Zoom in with Ctrl + P
		InputMap imap = table.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
		ActionMap amap = table.getActionMap();
		imap.put(KeyStroke.getKeyStroke("ctrl released P"), "zoomIn");
		amap.put("zoomIn", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				zoomChange(true);
			}

		});
		//Zoom in with Ctrl + O
		imap.put(KeyStroke.getKeyStroke("ctrl released M"), "zoomOut");
		amap.put("zoomOut", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				zoomChange(false);
			}

		});

		updateTable();
		updateViewState();

		// int WIDTH = getPreferredSize().width;
		// double scale = tableModel.getInitialZoomScale(
		// //jspContent.getVisibleRect().width,
		// WIDTH,
		// table.getIntercellSpacing().width);
		// setZoom((int) (scale * 100));
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
		DEBUG.printMessage(" *** Updating table: " + update_counter);// System.currentTimeMillis()
																		// );
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
			table.revalidate();
			// repaint();
		}

		jsZoom.setValue(data.zoom);
		isUpdating = updatingState;
		// }
	}// function updateTable

	private void startUpdater() {

		if (tableUpdaterTask != null)
			tableUpdaterTask.cancel();
		tableUpdaterTask = new TableUpdater();
		updateTimer.schedule(tableUpdaterTask, DELAY_TO_UPDATE);
	}// function startUpdater

	/**
	 * Adjust the columns depending on the current windows width and cell width
	 */
	private void adjustColumns() {
		int w = getSize().width;
		gallery.adjustColumn(w - 50);
		// DEBUG.printMessage(String.format(
		// "==>> JFrameGallery.autoAdjust: width: %d", w));
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
				e.printStackTrace();
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
		XmippImageWindow winAvg = new XmippImageWindow(impAvg, "AVG: "
				+ data.filename);
		WindowUtil.setLocation(0.2f, 0.5f, winAvg, this);
		winAvg.setVisible(true);
		XmippImageWindow winStd = new XmippImageWindow(impStd, "STD: "
				+ data.filename);
		WindowUtil.setLocation(0.8f, 0.5f, winStd, this);
		winStd.setVisible(true);
	}

	private void openAsStack() {
		// ImagesWindowFactory.openGalleryAsImagePlus(gallery);
	}

	private void openAs3D() {
		// ImagesWindowFactory.openGalleryAs3D(gallery);
	}

	private void save() {
		if (dlgSave == null) {
			dlgSave = new SaveJDialog(this);
			saveAs();
		} else {
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

	// Sets path and filename automatically.
	// String filename = gallery.getFilename() != null ? gallery
	// .getFilename() : "";
	// fc.setSelectedFile(new File(forceExtension(filename, ".xmd")));
	//
	// if (fc.showSaveDialog(this) != JFileChooser.CANCEL_OPTION) {
	// boolean response = true;
	// if (fc.getSelectedFile().exists()) {
	// response = JOptionPane.showConfirmDialog(null,
	// Labels.MESSAGE_OVERWRITE_FILE,
	// Labels.MESSAGE_OVERWRITE_FILE_TITLE,
	// JOptionPane.OK_CANCEL_OPTION,
	// JOptionPane.QUESTION_MESSAGE) == JOptionPane.OK_OPTION;
	// }
	//
	// if (response) {
	// String path = fc.getSelectedFile().getAbsolutePath();
	// if (gallery.saveAsMetadata(path, all)) {
	// JOptionPane.showMessageDialog(this,
	// Labels.MESSAGE_FILE_SAVED + path,
	// Labels.MESSAGE_FILE_SAVED_TITLE,
	// JOptionPane.INFORMATION_MESSAGE);
	// }
	// }
	// }

	private void saveAs() {
		dlgSave.setMdFilename(data.getMdFilename());
		if (dlgSave.showDialog()) {
			try {
				String path = dlgSave.getMdFilename();
				DEBUG.printMessage(path);
				if (dlgSave.isAppendMode())
					data.md.writeBlock(path);
				else
					data.md.write(path);
				if (dlgSave.doSaveImages()) {
					data.md.writeImages(dlgSave.getOutput(),
							dlgSave.isOutputIndependent(),
							dlgSave.getImageLabel());
				}
			} catch (Exception e) {
				XmippDialog.showError(this, e.getMessage());
			}
		}
		// Sets path and filename automatically.
		// String filename = gallery.getFilename() != null ? gallery
		// .getFilename() : "";
		// fc.setSelectedFile(new File(forceExtension(filename, ".stk")));
		//
		// if (fc.showSaveDialog(this) != JFileChooser.CANCEL_OPTION) {
		// boolean response = true;
		// if (fc.getSelectedFile().exists()) {
		// response = JOptionPane.showConfirmDialog(null,
		// Labels.MESSAGE_OVERWRITE_FILE,
		// Labels.MESSAGE_OVERWRITE_FILE_TITLE,
		// JOptionPane.OK_CANCEL_OPTION,
		// JOptionPane.QUESTION_MESSAGE) == JOptionPane.OK_OPTION;
		// }
		//
		// if (response) {
		// String path = fc.getSelectedFile().getAbsolutePath();
		// if (gallery.saveAsStack(path, all)) {
		// JOptionPane.showMessageDialog(this,
		// Labels.MESSAGE_FILE_SAVED + path,
		// Labels.MESSAGE_FILE_SAVED_TITLE,
		// JOptionPane.INFORMATION_MESSAGE);
		// }
		// }
		// }
	}

	static String forceExtension(String filename, String ext) {
		int dot = filename.lastIndexOf(".");
		return filename.substring(0, dot) + ext;
	}

	public void pca() {
		try {
			ImageGeneric image = new ImageGeneric();
			data.md.getPCAbasis(image);
			ImagePlus imp = XmippImageConverter.convertToImagePlus(image);
			// new XmippImageWindow(imp, "PCA: " + data.filename);
			imp.setTitle("PCA: " + data.filename);
			ImagesWindowFactory.captureFrame(imp);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public void fsc() {
		try {
			JFrameFSC frame = new JFrameFSC(data);
			WindowUtil.centerWindows(frame, this);
			frame.setVisible(true);
		} catch (Exception ex) {
			DEBUG.printException(ex);
		}
	}

	private void reslice(RESLICE_MODE mode) throws Exception {
		String command = null;

		switch (mode) {
		case TOP_Y:
			command = "Top";
			break;
		case RIGHT_X:
			command = "Right";
			break;
		}

		// Get volume ImagePlus.
		// String filename = gallery.getFilename();
		// ImagePlus volume = XmippImageConverter.loadImage(filename);
		//
		// // Reslice.
		// IJ.run(volume, "Reslice [/]...", "slice=1.000 start=" + command);
		// volume = WindowManager.getCurrentImage();
		// volume.getWindow().setVisible(false);
		//
		// // Save temp file.
		// int index = filename.lastIndexOf(".");
		// String name = filename.substring(
		// filename.lastIndexOf(File.separator) + 1, index);
		// String ext = filename.substring(index);
		// File f = File.createTempFile(name + "_" + command, ext);
		// f.deleteOnExit();
		//
		// XmippImageConverter.saveImage(volume, f.getAbsolutePath());
		// volume.close();
		//
		// // Open as gallery.
		// ImagesWindowFactory.openFileAsGallery(f.getCanonicalPath());
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
			boolean hasRender = data.hasRenderLabel();
			btnChangeView.setEnabled(hasRender);
			jsZoom.setEnabled(hasRender);
			jlZoom.setEnabled(hasRender);
		}
	}

	/** Reload table data */
	protected void reloadTableData() {
		try {
			table.removeAll();
			createModel();
			// gallery.setShowLabels(menu.getShowLabel());
			createTable();
			adjustColumns();
			menu.update();
			updateCombos();
			if (dlgSave != null)
				dlgSave.setInitialValues();
			setTitle(gallery.getTitle());
		} catch (Exception e) {
			e.printStackTrace();
		}
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
						jcbAutoAdjustColumnsStateChanged(evt);
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
				return data.selectedVol;
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
		// Dimension dim = evt.getComponent().getSize();
		// DEBUG.printMessage(dim.toString());
		// DEBUG.printMessage(evt.getComponent().getName());
		if (!isUpdating && autoAdjustColumns)
			adjustColumns();
	}

	private void tableMouseClicked(MouseEvent evt) {
		final Point p = evt.getPoint();
		int view_row = table.rowAtPoint(p);
		int view_col = table.columnAtPoint(p);

		if (evt.getButton() == MouseEvent.BUTTON1) { // Left click.
			if (evt.getClickCount() > 1) {
				Object item = table.getValueAt(view_row, view_col);

				if (item instanceof ImageItem) {
					try {
						DEBUG.printMessage(Filename.getXmippPath());
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					new XmippImageWindow(((ImageItem) item).getImage());
					ImagesWindowFactory.captureFrame(((ImageItem) item)
							.getImage());
				}
			} else {
				// Ctrl adds items to selection, otherwise previous ones are
				// removed.
				if (!evt.isControlDown()) {
					gallery.clearSelection();
				}

				if (evt.isShiftDown()) {
					gallery.touchRange(previousSelectedRow,
							previousSelectedCol, view_row, view_col);
				} else {
					gallery.touchItem(view_row, view_col);
				}
				isUpdating = true;
				jsGoToImage.setValue(gallery.getIndex(view_row, view_col) + 1);
				isUpdating = false;
				// table.repaint();
			}

			if (!evt.isShiftDown()) {
				previousSelectedRow = view_row;
				previousSelectedCol = view_col;
			}
		} else if (evt.getButton() == MouseEvent.BUTTON3) { // Right click.
			if (gallery.getSelectedCount() < 2) {
				gallery.clearSelection();
				gallery.touchItem(view_row, view_col);

				table.setRowSelectionInterval(view_row, view_row);
				table.setColumnSelectionInterval(view_col, view_col);

				// table.repaint();
			}

			final MouseEvent me = evt;
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					jpopUpMenuTable.show(me.getComponent(), p);
				}
			});
		}
	}

	private void jcbAutoAdjustColumnsStateChanged(ActionEvent evt) {// GEN-FIRST:event_jcbAutoAdjustColumnsStateChanged
		setAutoAdjustColumns(jcbAutoAdjustColumns.isSelected());
		if (autoAdjustColumns) {
			adjustColumns();
		}
	}

	class TableUpdater extends TimerTask {

		@Override
		public void run() {
			DEBUG.printMessage("Starting table update...");
			try {
				Thread.sleep(5000);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			updateTable();
			tableUpdaterTask = null;
		}
	}

	class GalleryMenu extends XmippMenuCreator {
		public final String FILE = "File";
		public final String FILE_OPEN = "File.Open_mi";
		public final String FILE_OPENWITH_IJ = "File.OpenWithIJ_mi";
		public final String FILE_OPENWITH_CHIMERA = "File.OpenWithChimera_mi";
		public final String FILE_SAVE = "File.Save_mi";
		public final String FILE_SAVEAS = "File.SaveAs_mi";
		public final String FILE_EXIT = "File.Exit_mi";
		public final String DISPLAY = "Display";
		public final String DISPLAY_NORMALIZE = "Display.Normalize_cb";
		public final String DISPLAY_SHOWLABELS = "Display.ShowLabels_cb";
		public final String DISPLAY_RENDERIMAGES = "Display.RenderImages_cb";
		public final String DISPLAY_APPLYGEO = "Display.ApplyGeo_cb";
		public final String DISPLAY_COLUMNS = "Display.Columns_mi";
		public final String DISPLAY_RESLICE = "Display.Reslice";
		public final String DISPLAY_RESLICE_X = "Display.Reslice.X_rb";
		public final String DISPLAY_RESLICE_Y = "Display.Reslice.Y_rb";
		public final String DISPLAY_RESLICE_Z = "Display.Reslice.Z_rb";
		public final String STATS = "Stats";
		public final String STATS_AVGSTD = "Stats.AvgStd_mi";
		public final String STATS_PCA = "Stats.Pca_mi";
		public final String STATS_FSC = "Stats.Fsc_mi";

		public GalleryMenu() {	
			super(MenuType.JMENUBAR);
			update();
		}// constructor GalleryMenu

		public void update() {
			boolean galMode = data.isGalleryMode();
			boolean volMode = data.isVolumeMode();
			setItemEnabled(FILE_OPENWITH_CHIMERA, volMode);
			setItemEnabled(FILE_SAVE, !volMode);
			setItemEnabled(FILE_SAVEAS, !volMode);
			setItemEnabled(DISPLAY_NORMALIZE, gallery.getNormalized());
			setItemEnabled(DISPLAY_APPLYGEO, data.containsGeometryInfo());
			setItemSelected(DISPLAY_APPLYGEO, data.useGeo);
			setItemEnabled(DISPLAY_RENDERIMAGES, !galMode && data.hasRenderLabel());
			setItemSelected(DISPLAY_RENDERIMAGES, data.globalRender);
			setItemEnabled(DISPLAY_COLUMNS, !galMode);
			setItemEnabled(DISPLAY_RESLICE, volMode);
			setItemEnabled(STATS, !volMode);
		}// function update

		
		@Override
		protected void createItems() throws Exception {
			//File
			addItem(FILE, "File");
			addItem(FILE_OPEN, "Open ...", null, "control released O");
			addItem(FILE_OPENWITH_IJ, "Open with ImageJ", "ij.gif", "control released J");
			addItem(FILE_OPENWITH_CHIMERA, "Open with Chimera", "chimera.gif", "control released H");
			addSeparator(FILE);
			addItem(FILE_SAVE, "Save", "save.gif", "control released S");
			addItem(FILE_SAVEAS, "Save as", "save_as.gif");
			addSeparator(FILE);
			addItem(FILE_EXIT, "Exit", null, "control released Q");
			//Display
			addItem(DISPLAY, "Display");
			addItem(DISPLAY_NORMALIZE, "Global normalization", null, "control released N");
			addItem(DISPLAY_SHOWLABELS, "Show labels", null, "control released L");
			addSeparator(DISPLAY);
			addItem(DISPLAY_RENDERIMAGES, "Render images", null, "control released R");
			addItem(DISPLAY_APPLYGEO, "Apply geometry", null, "control released G");
			addItem(DISPLAY_COLUMNS, "Columns ...", "columns.gif");
			addItem(DISPLAY_RESLICE, "Reslice");
			addItem(DISPLAY_RESLICE_X, "X");
			addItem(DISPLAY_RESLICE_Y, "Y");
			addItem(DISPLAY_RESLICE_Z, "Z");
			//Statistics
			addItem(STATS, "Statistics");
			addItem(STATS_AVGSTD, "Avg & Std images");
			addItem(STATS_PCA, "PCA");
			addItem(STATS_FSC, "FSC");			
		}//function createItems
		
		@Override
		protected void handleActionPerformed(ActionEvent e) {
			String cmd = e.getActionCommand();
			
			if (cmd.equals(DISPLAY_NORMALIZE)) {
				gallery.setNormalized(getItemSelected(DISPLAY_NORMALIZE));
			} else if (cmd.equals(DISPLAY_APPLYGEO)) {
				if (data.containsGeometryInfo()) {
					((MetadataGallery) gallery).setUseGeometry(getItemSelected(DISPLAY_APPLYGEO));
				}
			} else if (cmd.equals(DISPLAY_SHOWLABELS)) {
				gallery.setShowLabels(getItemSelected(DISPLAY_SHOWLABELS));
			} else if (cmd.equals(DISPLAY_RENDERIMAGES)) {
				gallery.setRenderImages(getItemSelected(DISPLAY_RENDERIMAGES));
			} else if (cmd.equals(DISPLAY_COLUMNS)) {
				ColumnsJDialog dialog = new ColumnsJDialog(JFrameGallery.this);
				boolean result = dialog.showDialog();
				if (result) {
					ArrayList<ColumnInfo> columns = dialog.getColumnsResult();
					isUpdating = true;
					((MetadataGallery) gallery).updateColumnInfo(columns);
					gallery.fireTableDataChanged();	
					setItemEnabled(DISPLAY_RENDERIMAGES, data.globalRender);
					//menu.enableRenderImages(data.globalRender);
					isUpdating = false;
				}
			} else if (cmd.equals(STATS_AVGSTD))
				runInBackground(Worker.STATS);
			else if (cmd.equals(STATS_PCA))
				runInBackground(Worker.PCA);
			else if (cmd.equals(STATS_FSC))
				runInBackground(Worker.FSC);
			else if (cmd.equals(FILE_OPEN)){
				if (fc.showOpenDialog(JFrameGallery.this) != JFileChooser.CANCEL_OPTION){
					if (fc.getSelectedFile().exists())
						ImagesWindowFactory.openFileAsDefault(fc.getSelectedFile().getPath());
					else 
						XmippDialog.showError(JFrameGallery.this, 
								String.format("File: '%s' doesn't exist.", fc.getSelectedFile().getPath()));
				}
				// String filename = gallery.getFilename() != null ? gallery
				// .getFilename() : "";
				// fc.setSelectedFile(new File(forceExtension(filename, ".stk")));
				//
				// if (fc.showSaveDialog(this) != JFileChooser.CANCEL_OPTION) {
				// boolean response = true;
				// if (fc.getSelectedFile().exists()) {
				// response = JOptionPane.showConfirmDialog(null,
				// Labels.MESSAGE_OVERWRITE_FILE,
				// Labels.MESSAGE_OVERWRITE_FILE_TITLE,
				// JOptionPane.OK_CANCEL_OPTION,
				// JOptionPane.QUESTION_MESSAGE) == JOptionPane.OK_OPTION;
				// }
				//
				// if (response) {
				// String path = fc.getSelectedFile().getAbsolutePath();
				// if (gallery.saveAsStack(path, all)) {
				// JOptionPane.showMessageDialog(this,
				// Labels.MESSAGE_FILE_SAVED + path,
				// Labels.MESSAGE_FILE_SAVED_TITLE,
				// JOptionPane.INFORMATION_MESSAGE);
				// }
			}
			else if (cmd.equals(FILE_SAVE)) {
				save();
			} else if (cmd.equals(FILE_SAVEAS)) {
				saveAs();
			} else if (cmd.equals(FILE_EXIT)) {
				System.exit(0);
			} else if (cmd.equals(FILE_OPENWITH_CHIMERA)) {
				try {
					String args = data.selectedVol;
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
			}
		}//function handleActionPerformed
	}// class GalleryMenu

	class GalleryPopupMenu extends XmippMenuCreator {
		protected Point location;
		
		public final static String ENABLED = "Enabled_mi";
		public final static String DISABLED = "Disabled_mi";
		public final static String SELECT = "Select";
		public final static String SELECT_ALL = "Select.All_mi";
		public final static String SELECT_TOHERE = "Select.ToHere_mi";
		public final static String SELECT_FROMHERE = "Select.FromHere_mi";


		@Override
		protected void createItems() throws Exception {
			addItem(ENABLED, "Enable");
			addItem(DISABLED, "Disable");
			addSeparator();
			addItem(SELECT, "Select");
			addItem(SELECT_ALL, "All", null, "control released A");
			addItem(SELECT_TOHERE, "To here");
			addItem(SELECT_FROMHERE, "From here");			
		}//function createItems
		
		public GalleryPopupMenu() {
			super(MenuType.JPOPUPMENU);
//			add(jmiEnable);
//			add(jmiDisable);
//			addSeparator();
//			add(jmiSelectFrom);
//			add(jmiSelectTo);
//
//			jmiEnable.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E,
//					InputEvent.CTRL_DOWN_MASK));
//			jmiDisable.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_D,
//					InputEvent.CTRL_DOWN_MASK));
			// jmiEnableFrom.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F,
			// InputEvent.SHIFT_DOWN_MASK));
			// jmiDisableFrom.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F,
			// InputEvent.CTRL_DOWN_MASK));
			// jmiEnableTo.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T,
			// InputEvent.SHIFT_DOWN_MASK));
			// jmiDisableTo.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T,
			// InputEvent.CTRL_DOWN_MASK));

			// jmiEnable.addActionListener(new ActionListener() {
			//
			// public void actionPerformed(ActionEvent e) {
			// gallery.enableSelectedItems();
			// }
			// });
			//
			// jmiDisable.addActionListener(new ActionListener() {
			//
			// public void actionPerformed(ActionEvent e) {
			// gallery.disableSelectedItems();
			// }
			// });
			//
			// jmiEnableAll.addActionListener(new ActionListener() {
			//
			// public void actionPerformed(ActionEvent e) {
			// gallery.enableAllItems();
			// }
			// });
			//
			// jmiDisableAll.addActionListener(new ActionListener() {
			//
			// public void actionPerformed(ActionEvent e) {
			// gallery.disableAllItems();
			// }
			// });
			//
			// jmiEnableFrom.addActionListener(new ActionListener() {
			//
			// public void actionPerformed(ActionEvent e) {
			// int row = table.rowAtPoint(location);
			// int col = table.columnAtPoint(location);
			//
			// gallery.setEnabledFrom(row, col, true);
			// }
			// });
			//
			// jmiEnableTo.addActionListener(new ActionListener() {
			//
			// public void actionPerformed(ActionEvent e) {
			// int row = table.rowAtPoint(location);
			// int col = table.columnAtPoint(location);
			//
			// gallery.setEnabledTo(row, col, true);
			// }
			// });
			//
			// jmiDisableFrom.addActionListener(new ActionListener() {
			//
			// public void actionPerformed(ActionEvent e) {
			// int row = table.rowAtPoint(location);
			// int col = table.columnAtPoint(location);
			//
			// gallery.setEnabledFrom(row, col, false);
			// }
			// });
			//
			// jmiDisableTo.addActionListener(new ActionListener() {
			//
			// public void actionPerformed(ActionEvent e) {
			// int row = table.rowAtPoint(location);
			// int col = table.columnAtPoint(location);
			//
			// gallery.setEnabledTo(row, col, false);
			// }
			// });
		}// constructor JPopUpMenuGallery

		public void show(Component cmpnt, Point location) {
			this.location = location;

			// Update menu items status depending on item.
			int row = table.rowAtPoint(location);
			int col = table.columnAtPoint(location);
			// boolean enabled = ((I) table.getValueAt(row,
			// col)).isEnabled();

			// jmiDisable.setEnabled(enabled);
			// jmiEnable.setEnabled(!enabled);

			getPopupMenu().show(cmpnt, location.x, location.y);
		}// function show

		@Override
		protected void handleActionPerformed(ActionEvent evt) {
			// TODO Auto-generated method stub
			
		}

	}// class JPopUpMenuGallery
}// class JFrameGallery
