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

import xmipp.viewer.windows.ImagesWindowFactory;

import xmipp.ij.XmippImageConverter;
import xmipp.ij.XmippImageWindow;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippLabel;
import xmipp.utils.Param;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageWindow;
import xmipp.utils.XmippResource;

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
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Timer;
import java.util.TimerTask;

import javax.swing.AbstractButton;
import javax.swing.ButtonGroup;
import javax.swing.ComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
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
	JFileChooser fc = new JFileChooser();

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

		jspContent = new javax.swing.JScrollPane();
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
		setJMenuBar(menu);
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

		updateTable();
		updateViewState();

		// int WIDTH = getPreferredSize().width;
		// double scale = tableModel.getInitialZoomScale(
		// //jspContent.getVisibleRect().width,
		// WIDTH,
		// table.getIntercellSpacing().width);
		// setZoom((int) (scale * 100));
	}// function createTable

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

	private void avgImage() {
		try {
			ImageGeneric imgAvg = new ImageGeneric();
			ImageGeneric imgStd = new ImageGeneric();
			data.md.getStatsImages(imgAvg, imgStd, true);
			ImagePlus impAvg = XmippImageConverter.convertImageGenericToImageJ(imgAvg);
			ImagePlus impStd =  XmippImageConverter.convertImageGenericToImageJ(imgStd);
			imgAvg.destroy();
			imgStd.destroy();
			impAvg.setTitle("AVG: " + data.filename);
			impStd.setTitle("STD: " + data.filename);
			 ImagesWindowFactory.captureFrame(impAvg);
			 ImagesWindowFactory.captureFrame(impStd);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// ImagesWindowFactory.captureFrame(ImageOperations.mean(tableModel.getAllItems()));
	}

	private void stdDevImage() {
		// ImagesWindowFactory.captureFrame(ImageOperations.std_deviation(tableModel.getAllItems()));
	}

	private void openAsStack() {
		// ImagesWindowFactory.openGalleryAsImagePlus(gallery);
	}

	private void openAs3D() {
		// ImagesWindowFactory.openGalleryAs3D(gallery);
	}

	private void save() {
		SaveJDialog dlg = new SaveJDialog(this);
		if (dlg.showDialog())
		{
			try {
				String path = dlg.getData();
				DEBUG.printMessage(path);
				data.md.write(path);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
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
	}

	private void saveAsStack(boolean all) {
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
		 ImagePlus imp = XmippImageConverter.convertImageGenericToImageJ(image);
		 imp.setTitle("PCA: " + data.filename);
		 ImagesWindowFactory.captureFrame(imp);
		 } catch (Exception ex) {
		 DEBUG.printException(ex);
		 }
	}

	public void fsc() {
		// String filename = gallery.getFilename();
		//
		 try {
		 ImagesWindowFactory.openFSCWindow(data.md);
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

		// TODO: MOVE THIS TO MENU OPTION
		// jcbShowLabels = new javax.swing.JCheckBox();
		// jcbMDLabels = new javax.swing.JComboBox();
		// jcbShowLabels.setText(Labels.LABEL_SHOW_LABELS);
		// jcbShowLabels.addActionListener(new java.awt.event.ActionListener() {
		// public void actionPerformed(java.awt.event.ActionEvent evt) {
		// jcbShowLabelsActionPerformed(evt);
		// }
		// });
		// toolBar.add(jcbShowLabels);
		//
		// jcbMDLabels.addItemListener(new java.awt.event.ItemListener() {
		// public void itemStateChanged(java.awt.event.ItemEvent evt) {
		// jcbMDLabelsItemStateChanged(evt);
		// }
		// });
		// toolBar.add(jcbMDLabels);
		//
		// jcbSortByLabel.setText(Labels.LABEL_SORT_BY_LABEL);
		// jcbSortByLabel.addActionListener(new java.awt.event.ActionListener()
		// {
		// public void actionPerformed(java.awt.event.ActionEvent evt) {
		// jcbSortByLabelActionPerformed(evt);
		// }
		// });
		// Set comboBox items.
		// String labels[] = tableModel.getLabels();
		// for (int i = 0; i < labels.length; i++) {
		// jcbMDLabels.addItem(labels[i]);
		// }
		// toolBar.add(jcbSortByLabel);

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

	private void jsGoToImageStateChanged(javax.swing.event.ChangeEvent evt) {// GEN-FIRST:event_jsGoToImageStateChanged
		DEBUG.printMessage("jsGotoImage...");
		if (!isUpdating) {
			goToImage((Integer) jsGoToImage.getValue() - 1);
		}
	}

	private void formComponentResized(java.awt.event.ComponentEvent evt) {// GEN-FIRST:event_formComponentResized
		// Dimension dim = evt.getComponent().getSize();
		// DEBUG.printMessage(dim.toString());
		// DEBUG.printMessage(evt.getComponent().getName());
		if (!isUpdating && autoAdjustColumns)
			adjustColumns();
	}

	private void tableMouseClicked(java.awt.event.MouseEvent evt) {// GEN-FIRST:event_tableMouseClicked
		final Point p = evt.getPoint();
		int view_row = table.rowAtPoint(p);
		int view_col = table.columnAtPoint(p);

		if (evt.getButton() == MouseEvent.BUTTON1) { // Left click.
			if (evt.getClickCount() > 1) {
				Object item = table.getValueAt(view_row, view_col);

				if (item instanceof ImageItem) {
//					ImagesWindowFactory.captureFrame(((ImageItem) item)
//							.getImage());
					new XmippImageWindow(((ImageItem) item).getImage());
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

	private void jtbNormalizeActionPerformed(java.awt.event.ActionEvent evt) {// GEN-FIRST:event_jtbNormalizeActionPerformed
		// setNormalized(jtbNormalize.isSelected());
		// updateTable();
	}// GEN-LAST:event_jtbNormalizeActionPerformed

	private void jcbMDLabelsItemStateChanged(java.awt.event.ItemEvent evt) {// GEN-FIRST:event_jcbMDLabelsItemStateChanged
		if (evt.getStateChange() == ItemEvent.DESELECTED) {
			// gallery.setSelectedLabel(jcbMDLabels.getSelectedIndex());
			updateTable();
		}
	}// GEN-LAST:event_jcbMDLabelsItemStateChanged

	private void jcbSortByLabelActionPerformed(java.awt.event.ActionEvent evt) {// GEN-FIRST:event_jcbSortByLabelActionPerformed
		// gallery.setSorting(jcbSortByLabel.isSelected());
		updateTable();
	}// GEN-LAST:event_jcbSortByLabelActionPerformed

	private void jtbUseGeometryActionPerformed(java.awt.event.ActionEvent evt) {// GEN-FIRST:event_jtbUseGeometryActionPerformed
		// if (tableModel.containsGeometryInfo()) {
		// ((MDTableModel)
		// tableModel).setUseGeometry(jtbUseGeometry.isSelected());
		// updateTable();
		// }
	}// GEN-LAST:event_jtbUseGeometryActionPerformed

	private void jcbAutoAdjustColumnsStateChanged(ActionEvent evt) {// GEN-FIRST:event_jcbAutoAdjustColumnsStateChanged
		setAutoAdjustColumns(jcbAutoAdjustColumns.isSelected());
		if (autoAdjustColumns) {
			adjustColumns();
		}
	}// GEN-LAST:event_jcbAutoAdjustColumnsStateChanged
		// Variables declaration - do not modify//GEN-BEGIN:variables

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

	class GalleryMenu extends JMenuBar implements ActionListener {
		protected JMenu jmFile = new JMenu(XmippLabel.LABEL_GALLERY_FILE);
		protected JMenu jmSave222 = new JMenu(XmippLabel.LABEL_GALLERY_SAVE);
		protected JMenu jmDisplay = new JMenu("Display");

		protected JMenuItem jmiSave = new JMenuItem(
				XmippLabel.LABEL_GALLERY_SAVE, XmippResource.getIcon("save.gif"));
		protected JMenuItem jmiSaveAs = new JMenuItem(
				XmippLabel.LABEL_GALLERY_SAVEAS, XmippResource.getIcon("save_as.gif"));
//		protected JMenuItem jmiSaveSelectionAsMetadata = new JMenuItem(
//				XmippLabel.LABEL_GALLERY_SAVE_SELECTION_AS_METADATA);
//		protected JMenuItem jmiSaveSelectionAsStack = new JMenuItem(
//				XmippLabel.LABEL_GALLERY_SAVE_SELECTION_AS_IMAGE);
		protected JMenuItem jmiExit = new JMenuItem(
				XmippLabel.LABEL_GALLERY_EXIT);
		protected JMenu jmStatistics = new JMenu(XmippLabel.MENU_STATS);
		protected JMenuItem jmiAVG = new JMenuItem(XmippLabel.BUTTON_MEAN);
		protected JMenuItem jmiSTDEV = new JMenuItem(
				XmippLabel.BUTTON_STD_DEVIATION);
		protected JMenuItem jmiPCA = new JMenuItem(XmippLabel.BUTTON_PCA);
		protected JMenuItem jmiFSC = new JMenuItem(XmippLabel.BUTTON_FSC);
		protected JMenu jmOpenWith = new JMenu(XmippLabel.MENU_OPEN_WITH);
		protected JMenuItem jmiOpenWithChimera = new JMenuItem(
				XmippLabel.MENU_OPEN_WITH_CHIMERA,
				XmippResource.getIcon("chimera.gif"));
		protected JMenuItem jmiOpenWithImageJ = new JMenuItem(
				XmippLabel.MENU_OPEN_WITH_IJ, XmippResource.getIcon("ij.gif"));
		// protected JMenu jmReslice = new JMenu(Labels.LABEL_RESLICE);
		// protected JMenuItem jmiResliceTop = new JMenuItem(
		// Labels.LABEL_RESLICE_TOP);
		// protected JMenuItem jmiResliceRight = new JMenuItem(
		// Labels.LABEL_RESLICE_RIGHT);

		protected JCheckBoxMenuItem jmiNormalize = new JCheckBoxMenuItem(
				XmippLabel.MSG_NORMALIZE);
		protected JCheckBoxMenuItem jmiApplyGeo = new JCheckBoxMenuItem(
				XmippLabel.MSG_APPLY_GEO);
		protected JMenuItem jmiColumns = new JMenuItem("Columns ...",
				XmippResource.getIcon("columns.gif"));
		protected JCheckBoxMenuItem jmiShowLabel = new JCheckBoxMenuItem(
				XmippLabel.MSG_SHOW_LABEL);
		protected JCheckBoxMenuItem jmiRenderImage = new JCheckBoxMenuItem(
				XmippLabel.MSG_RENDER_IMG);
		protected JMenu jmReslice = new JMenu(XmippLabel.LABEL_RESLICE);
		protected JRadioButtonMenuItem jmiAxisX = new JRadioButtonMenuItem("X");
		protected JRadioButtonMenuItem jmiAxisY = new JRadioButtonMenuItem("Y");
		protected JRadioButtonMenuItem jmiAxisZ = new JRadioButtonMenuItem("Z");

		// public void enableNormalize(boolean value) {
		// jmiNormalize.setEnabled(value);
		// }
		//
		// public void selectNormalize(boolean value) {
		// jmiNormalize.setSelected(value);
		// }

		public boolean getNormalize() {
			return jmiNormalize.isSelected();
		}

		// public void enableApplyGeo(boolean value) {
		// jmiApplyGeo.setEnabled(value);
		// }
		//
		// public void selectApplyGeo(boolean value) {
		// jmiApplyGeo.setSelected(value);
		// }

		public boolean getApplyGeo() {
			return jmiApplyGeo.isSelected();
		}

		// public void enableShowLabel(boolean value) {
		// jmiShowLabel.setEnabled(value);
		// }

		public boolean getShowLabel() {
			return jmiShowLabel.isSelected();
		}

		 public void enableRenderImages(boolean value) {
		 jmiRenderImage.setEnabled(value);
		 }

		public boolean getRenderImages() {
			return jmiRenderImage.isSelected();
		}

		public void enableReslice(boolean value) {
			jmReslice.setEnabled(value);
		}

		public void actionPerformed(ActionEvent e) {
			JMenuItem jmi = (JMenuItem) e.getSource();
			if (jmi == jmiNormalize) {
				gallery.setNormalized(jmiNormalize.isSelected());
			} else if (jmi == jmiApplyGeo) {
				if (data.containsGeometryInfo()) {
					((MetadataGallery) gallery).setUseGeometry(jmiApplyGeo
							.isSelected());
				}
			} else if (jmi == jmiShowLabel) {
				gallery.setShowLabels(jmiShowLabel.isSelected());
			} else if (jmi == jmiRenderImage) {
				gallery.setRenderImages(jmiRenderImage.isSelected());
			} else if (jmi == jmiColumns) {
				ColumnsJDialog dialog = new ColumnsJDialog(JFrameGallery.this);
				DEBUG.printMessage("Before showing dialog");
				boolean result = dialog.showDialog();
				DEBUG.printMessage("result:" + (result ? "True": "False"));
				if (result) {
					DEBUG.printMessage("AFter showing dialog");
					ArrayList<ColumnInfo> columns = dialog.getColumnsResult();
					isUpdating = true;
					((MetadataGallery) gallery).updateColumnInfo(columns);
					gallery.fireTableDataChanged();
					menu.enableRenderImages(data.globalRender);
					isUpdating = false;
				}
			} else if (jmi == jmiAVG)
				avgImage();
			else if (jmi == jmiSTDEV)
				stdDevImage();
			else if (jmi == jmiPCA)
				pca();
			else if (jmi == jmiFSC)
				fsc();
			else if (jmi == jmiSave) {
				save();
//			} else if (jmi == jmiSaveSelectionAsMetadata) {
//				saveAsMetadata(false);
			} else if (jmi == jmiSaveAs) {
				saveAsStack(true);
//			} else if (jmi == jmiSaveSelectionAsStack) {
//				saveAsStack(false);
			} else if (jmi == jmiExit) {
				System.exit(0);
			} else if (jmi == jmiOpenWithChimera) {
				try {
					String args = data.selectedVol;
					if (Filename.isSpiderVolume(args))
						args = "spider:" + args;
					// FIXME: Check chimera is installed and volume is on spider
					// format
					Process p = new ProcessBuilder("chimera", args).start();
				} catch (Exception ex) {
					ex.printStackTrace();
				}
			} else if (jmi == jmiOpenWithImageJ) {
				try {
					ImagePlus imp = gallery.getImagePlus();
					ImagesWindowFactory.captureFrame(imp);
				} catch (Exception e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
			}
		}

		private void addMenuItem(JMenu parent, JMenuItem item){
			parent.add(item);
			item.addActionListener(this);
		}
		
		public GalleryMenu() {
			super();

			// File menu
//			jmSave222.addSeparator();
//			jmSave222.add(jmiSaveSelectionAsMetadata);
//			jmSave222.add(jmiSaveSelectionAsStack);

			add(jmFile);
			addMenuItem(jmFile, jmiSave);
			addMenuItem(jmFile, jmiSaveAs);
			//addMenuItem(jmFile, jmSave222);
			jmFile.addSeparator();
			addMenuItem(jmFile, jmiExit);

			// Display menu
			add(jmDisplay);
			addMenuItem(jmDisplay, jmiNormalize);
			addMenuItem(jmDisplay, jmiShowLabel);
			jmDisplay.addSeparator();
			addMenuItem(jmDisplay, jmiRenderImage);
			addMenuItem(jmDisplay, jmiApplyGeo);
			addMenuItem(jmDisplay, jmiColumns);
			jmDisplay.addSeparator();
			addMenuItem(jmDisplay, jmReslice);
			addMenuItem(jmReslice, jmiAxisX);
			addMenuItem(jmReslice, jmiAxisY);
			addMenuItem(jmReslice, jmiAxisZ);
			ButtonGroup group = new ButtonGroup();
			group.add(jmiAxisX);
			group.add(jmiAxisY);
			group.add(jmiAxisZ);

//			jmiNormalize.addActionListener(this);
//			jmiApplyGeo.addActionListener(this);
//			jmiShowLabel.addActionListener(this);
//			jmiRenderImage.addActionListener(this);
//			jmiColumns.addActionListener(this);

			// Statistics menu
			add(jmStatistics);
			addMenuItem(jmStatistics, jmiAVG);
			addMenuItem(jmStatistics, jmiSTDEV);
			addMenuItem(jmStatistics, jmiPCA);
			addMenuItem(jmStatistics, jmiFSC);

			// Open with menu
			add(jmOpenWith);
			addMenuItem(jmOpenWith, jmiOpenWithChimera);
			addMenuItem(jmOpenWith, jmiOpenWithImageJ);
			// addMenuItem(jmOpenWith, jmiOpenAsStack);

//			jmiOpenWithChimera.addActionListener(this);
//			jmiOpenWithImageJ.addActionListener(this);

			update();
		}

		public void update() {
			boolean galMode = data.isGalleryMode();
			boolean volMode = data.isVolumeMode();
			jmStatistics.setEnabled(!volMode);
			jmReslice.setEnabled(volMode);
			// jmTopAxis.setEnabled(isVolume);
			jmiOpenWithChimera.setEnabled(volMode);
			// jmiOpenAs3D.setEnabled(!isMetaData);
			// jmiOpenAsMetadata.setEnabled(!isVolume);

			// Volumes can't be saved as metadata.
			jmiSave.setEnabled(!volMode);
			jmiSaveAs.setEnabled(!volMode);
			//jmiSaveSelectionAsMetadata.setEnabled(!volMode);
			jmiApplyGeo.setEnabled(data.containsGeometryInfo());
			jmiApplyGeo.setSelected(data.useGeo);
			jmiNormalize.setEnabled(gallery.getNormalized());
			jmiRenderImage.setEnabled(!galMode && data.hasRenderLabel());
			jmiRenderImage.setSelected(data.globalRender);
			jmiColumns.setEnabled(!galMode);
		}
	}

	class GalleryPopupMenu extends JPopupMenu {

		protected Point location;
		protected JMenuItem jmiEnable = new JMenuItem(
				XmippLabel.LABEL_GALLERY_ENABLE);
		protected JMenuItem jmiDisable = new JMenuItem(
				XmippLabel.LABEL_GALLERY_DISABLE);
		protected JMenuItem jmiEnableAll = new JMenuItem(
				XmippLabel.LABEL_GALLERY_ENABLE_ALL);
		protected JMenuItem jmiDisableAll = new JMenuItem(
				XmippLabel.LABEL_GALLERY_DISABLE_ALL);
		protected JMenuItem jmiEnableFrom = new JMenuItem(
				XmippLabel.LABEL_GALLERY_ENABLE_FROM);
		protected JMenuItem jmiEnableTo = new JMenuItem(
				XmippLabel.LABEL_GALLERY_ENABLE_TO);
		protected JMenuItem jmiDisableFrom = new JMenuItem(
				XmippLabel.LABEL_GALLERY_DISABLE_FROM);
		protected JMenuItem jmiDisableTo = new JMenuItem(
				XmippLabel.LABEL_GALLERY_DISABLE_TO);

		public GalleryPopupMenu() {
			add(jmiEnable);
			add(jmiDisable);
			addSeparator();
			add(jmiEnableAll);
			add(jmiDisableAll);
			addSeparator();
			add(jmiEnableFrom);
			add(jmiEnableTo);
			addSeparator();
			add(jmiDisableFrom);
			add(jmiDisableTo);

			jmiEnable.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E,
					InputEvent.CTRL_DOWN_MASK));
			jmiDisable.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_D,
					InputEvent.CTRL_DOWN_MASK));
			jmiEnableAll.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_A,
					InputEvent.CTRL_DOWN_MASK));
			jmiDisableAll.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N,
					InputEvent.CTRL_DOWN_MASK));

			// jmiEnableFrom.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F,
			// InputEvent.SHIFT_DOWN_MASK));
			// jmiDisableFrom.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F,
			// InputEvent.CTRL_DOWN_MASK));
			// jmiEnableTo.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T,
			// InputEvent.SHIFT_DOWN_MASK));
			// jmiDisableTo.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T,
			// InputEvent.CTRL_DOWN_MASK));

			table.addKeyListener(new KeyListener() {

				public void keyTyped(KeyEvent e) {
				}

				public void keyPressed(KeyEvent e) {
				}

				public void keyReleased(KeyEvent e) {
					// if (e.isControlDown()) {
					// if (e.getKeyCode() == KeyEvent.VK_E) {
					// gallery.enableSelectedItems();
					// }
					// if (e.getKeyCode() == KeyEvent.VK_D) {
					// gallery.disableSelectedItems();
					// }
					// if (e.getKeyCode() == KeyEvent.VK_A) {
					// gallery.enableAllItems();
					// }
					// if (e.getKeyCode() == KeyEvent.VK_N) {
					// gallery.disableAllItems();
					// }
					// }
				}
			});

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

			show(cmpnt, location.x, location.y);
		}// function show
	}// class JPopUpMenuGallery

}// class JFrameGallery
