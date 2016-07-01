/***************************************************************************
 * Authors:     Juanjo Vega
 * 		J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.AbstractAction;
import javax.swing.AbstractButton;
import javax.swing.ActionMap;
import javax.swing.ButtonGroup;
import javax.swing.ComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
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
import javax.swing.SwingUtilities;
import javax.swing.event.ListDataListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.JTableHeader;

import xmipp.ij.commons.ImagePlusLoader;
import xmipp.ij.commons.Tool;
import xmipp.ij.commons.XmippApplication;
import xmipp.ij.commons.XmippUtil;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.Params;
import xmipp.utils.QuickHelpJDialog;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippLabel;
import xmipp.utils.XmippMenuBarCreator;
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
import xmipp.viewer.models.MetadataTableModel;
import xmipp.viewer.particlepicker.extract.ExtractParticlePicker;
import xmipp.viewer.particlepicker.extract.ExtractPickerJFrame;

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
	protected ImageGalleryTableModel gallery;
	private GalleryRowHeaderModel rowHeaderModel;
	private int previousSelectedRow, previousSelectedCol;
	private JList rowHeader;
	// this flag will be used to avoid firing properties change events
	// when the change is from our code and not external user interaction
	private boolean isUpdating;
	private GalleryPopupMenu jpopUpMenuTable;
	protected GalleryMenu menu;
	protected XmippFileChooser fc;
	private SaveJDialog dlgSave = null;
	private boolean saved = false;
	private ClassesJDialog dlgClasses = null;

	private JLabel jlZoom;
	private JLabel jlGoToImage;
	private JLabel jlRows;
	private JLabel jlColumns;
	private JToggleButton jcbAutoAdjustColumns;
	private JButton btnChangeView;
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
	protected static final int MIN_WIDTH = 800;
	protected static int MIN_HEIGHT;
	protected static int MAX_HEIGHT;
	protected static int MAX_WIDTH;
	protected static Dimension screenSize;
	


	/** Store data about visualization */
	protected GalleryData data;
	private ExtractPickerJFrame extractframe;
	private ButtonGroup reslicegroup;
	protected JPanel buttonspn;
	protected JButton searchbt;
    protected JButton plotbt;
	protected JButton chimerabt;

	/** Some static initialization for fancy default dimensions */
	static
	{
		screenSize = XmippWindowUtil.getScreenRectangle().getSize();
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
            data.setWindow(this);
			createModel();
			createGUI();
			XmippApplication.addInstance(false);
		}
		catch (Exception e)
		{
			DEBUG.printException(e);
            setVisible(false);
            dispose();
            
            throw new IllegalArgumentException(e);
		}
	}

	/** Constructors */
	public GalleryJFrame(String filename, Params parameters)
	{
		super();
		init(new GalleryData(this, filename, parameters));
	}

	public GalleryJFrame(MetaData md, Params parameters)
	{
		super();
		init(new GalleryData(this, parameters, md));
	}
        
        public GalleryJFrame(GalleryData data)
	{
		super();
		init(data);
	}

	/**
	 * Open another metadata separately *
	 */
	public void openMetadata(final MetaData md)
	{
            try
            {
                SwingUtilities.invokeLater(new Runnable() {

                    @Override
                    public void run() {
                        new GalleryJFrame(md, new Params());
                    }
                });
                
            }
            catch(Exception e)
            {
                XmippDialog.showError(this, e.getMessage());
            }

	}

	protected void createModel() throws Exception
	{
		createModel(null);
	}
	/**
	 * Function to create the gallery type depending on the filename
	 * 
	 * @throws Exception
	 */
	protected void createModel(boolean[] selection) throws Exception
	{
		gallery = data.createModel(selection);
        if (data.getModelColumns() != null)
			gallery.setColumns(data.getModelColumns());
		else if (data.getModelRows() != null)
			gallery.setRows(data.getModelRows());
        int index = gallery.getSelTo();
        if(jsGoToImage != null)
        	jsGoToImage.setValue(index + 1);
	}

	public GalleryData getData()
	{
		return data;
	}

	/** Close the application, check if changes first */
	public void close()
	{
		close(true);
	}// function close
        
        /** Close the application, check if changes first */
	public void close(boolean ask)
	{
            boolean isclose = true;
            if(ask)
                isclose = proceedWithChanges();
            if (isclose)
            {
                    setVisible(false);
                    dispose();
                    XmippApplication.removeInstance(false);

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
	protected void createGUI()
	{
		// Create file chooser and set current dir
		//setIconImage(XmippResource.getIcon("xmipp_logo.png").getImage());
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
		container.add(toolBar, XmippWindowUtil.getConstraints(c, 0, 0, 1, 1, GridBagConstraints.HORIZONTAL));
		setInitialValues();

		// Create combos for selection of blocks and/or volumes
		createCombos();
		container.add(cbPanel, XmippWindowUtil.getConstraints(c, 0, 1, 1, 1, GridBagConstraints.HORIZONTAL));
		updateVisibleCombos();

		jspContent = new GalleryScroll();
		// Create table
		createTable();
		
		c.weightx = 1.0;
		c.weighty = 1.0;
		container.add(jspContent, XmippWindowUtil.getConstraints(c, 0, 2, 1, 1, GridBagConstraints.BOTH));
		c.weightx = 0;
		c.weighty = 0;
		
		buttonspn = new JPanel(new FlowLayout(FlowLayout.TRAILING));
		container.add(buttonspn, XmippWindowUtil.getConstraints(c, 0, 3, 1, 1, GridBagConstraints.HORIZONTAL));
						
		// Create the menu for table
		initGalleryMenu();
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
				changeView();
			}
		});
		enableToolBarActions();
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
		XmippWindowUtil.setScipionImageIcon(this);
	}

	private void setInitialValues()
	{
		boolean adjust = true;
		if (data.parameters.columns > 0 || data.parameters.rows > 0)
			adjust = false;


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

        public void fireTableRowsUpdated(int from, int to) {
            gallery.fireTableRowsUpdated(from, to);
        }

    

        protected void initGalleryMenu() {
            menu = new GalleryMenu();
                    
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
			rowHeaderModel = (data.isColumnFormat() || !data.isTableMode()) ? new GalleryRowHeaderModel(gallery.getRowCount(), 1)
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
			public void mouseReleased(java.awt.event.MouseEvent evt)
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
                                if(table.getSelectedRow() == -1)
                                    return;
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

                case KeyEvent.VK_SPACE:
                        int from = gallery.getSelFrom(), to = gallery.getSelTo();
                        for (int i = from; i <= to; ++i)
                            if (gallery.isSelected(i))
                                data.setEnabled(i, !data.isEnabled(i));
                        if(from != -1)
                            gallery.fireTableRowsUpdated(from, to);
				}
				if (vdir != 0 || hdir != 0)
				{
					
					int newRow = table.getSelectedRow() + vdir;
					int col = (table.getSelectedColumn() == -1)? 0 : table.getSelectedColumn() + hdir;//col is -1 in metadata mode
					//System.out.printf("newRow: %d, col: %d\n", newRow, col);
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

		jsZoom.setValue(data.getZoom());
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
		if (data.isAutoAdjust()){
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
		jcbAutoAdjustColumns.setSelected(autoAdjustColumns);
		jcbAutoAdjustColumns.setEnabled(!data.isTableMode());
		jsColumns.setEnabled(!autoAdjustColumns);
		jsRows.setEnabled(!autoAdjustColumns);
	}

	private void makeVisible(int row, int col)
	{
		//System.out.println(String.format("gotoImage,  row: %d, col:%d", row, col));

		// Gets current selected cell bounds.
		Rectangle rect = table.getCellRect(row, col, true);
		// Ensures item is visible
		Point pos = jspContent.getViewport().getViewPosition();
		rect.translate(-pos.x, -pos.y);
		jspContent.getViewport().scrollRectToVisible(rect);

		repaint();
	}

	private void goToImage(int index)
	{
		gallery.gotoItem(index);
		Point coords = gallery.getCoords(index);
		makeVisible(coords.y, coords.x);
	}

	

	/** Function to create and launch the worker, blocking the gui */
	public void runInBackground(int operation)
	{

                boolean selected = false;
                if(gallery.hasSelection())
                {
                    int result = XmippDialog.showQuestionYesNoCancel(this, "This operation processes all images by default. Would you like to use selection instead?");
                    
                    if(result == XmippQuestionDialog.YES_OPTION)
                        selected = true;
                    else if(result == XmippQuestionDialog.CANCEL_OPTION)
                       
                        return;
                }
                
                
                
		Worker w = new Worker(operation, selected, this);
		XmippWindowUtil.blockGUI(this, w.getMessage());
		Thread thr = new Thread(w);
		thr.start();
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
		if (data.isTableMode())
		{ // if we are in table mode only allow change
			// if exist render label
			enableToolBarActions();
            autoAdjustColumns(true);
		}
        else
            autoAdjustColumns(data.isAutoAdjust());
	}
	public void reloadTableData()
	{
		reloadTableData(true, null);
	}

	public void reloadTableData(boolean[] selection)
	{
		reloadTableData(true, selection);
	}

	/** Reload table data */
	public void reloadTableData(boolean changed, boolean[] selection)
	{
		try
		{
			//DEBUG.printMessage("reloadTableData...");
			if (table != null)
				table.removeAll();
			createModel(selection);
			// gallery.setShowLabels(menu.getShowLabel());
			createTable();

			menu.update();
			updateVisibleCombos();
            searchbt.setEnabled(data.isTableMode());
            plotbt.setEnabled(data.isTableMode());
			if (dlgSave != null && changed)
				dlgSave.setInitialValues();

			this.saved = !changed;
			setGalleryTitle();
            Integer rows = data.getModelRows();
            Integer cols = data.getModelColumns();
			if (rows != null)
				gallery.setRows(rows);
			if (cols != null)
				gallery.setColumns(cols);
			

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
		reloadTableData(changed, gallery.getSelection());
		data.setMdChanges(changed);

	}// function reloadMd

	/**
	 * Fill some label mode can be: "constant", "linear", "uniform", "gaussian"
	 * values is a list of string
	 * */
	public void fillLabel(int label, String mode, String... values) throws Exception
	{
		if (mode.equalsIgnoreCase(MetaData.FILL_CONSTANT))
			data.fillConstant(label, values[0]);
		else
		{
			Double v1 = Double.parseDouble(values[0]);
			Double v2 = Double.parseDouble(values[1]);
			if (mode.equalsIgnoreCase(MetaData.FILL_LINEAR))
				data.fillLinear(label, v1, v2);
			else if (mode.equalsIgnoreCase(MetaData.FILL_RAND_UNIFORM))
				data.fillRandom(label, "uniform", v1, v2);
			else if (mode.equalsIgnoreCase(MetaData.FILL_RAND_GAUSSIAN))
				data.fillRandom(label, "gaussian", v1, v2);
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
					data.removeSelection(gallery.getSelection());
				else
					data.removeDisabled();
				gallery.clearSelection();
				reloadMd();
			}
		}
		catch (Exception ex)
		{
			showException(ex);
		}
	}

	

	/** Find and replace in metadata */
	public void findReplace() 
	{
		MDSearchJDialog dlg = new MDSearchJDialog(this, table, data);
		dlg.setVisible(true);
	}

	/** Drop some label from the metadata */
	public void removeLabel(int label) throws Exception
	{
		data.removeLabel(label);
		reloadMd();
	}

	/***
	 * Function to create the main toolbar
	 */
        
        protected void changeView()
        {
            	data.changeMode();
            	boolean[] selection = null;
            	if(!data.isVolumeMd())
            		selection = gallery.getSelection();
            	reloadTableData(selection);
            	if(gallery.getSelTo() != -1)
            	{
	            	Point coords = gallery.getCoords(gallery.getSelTo());
					makeVisible(coords.y, coords.x);
            	}
	
        }
        
	protected void createToolbar()
	{
		// Create Main TOOLBAR
		toolBar = new JToolBar();
                toolBar.setFloatable(false);
		toolBar.setRollover(true);
		toolBar.setLayout(new FlowLayout(FlowLayout.LEFT));

		btnChangeView = new JButton();
		// updateViewState();
		btnChangeView.addActionListener(new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				changeView();
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
				if (gallery.getCellSize().getHeight() < 30 && zoom < gallery.data.getZoom())
				{
					jsZoom.setValue(gallery.data.getZoom());//keep previous zoom
					return;
				}
				gallery.setZoom(zoom);
				makeVisible(gallery.getFirstSelectedIndex(), 0);
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
		jcbAutoAdjustColumns.setSelected(data.isAutoAdjust());
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

//		if (gallery.getSize() > 0)
//		{
//
//			jsRows.setModel(new SpinnerNumberModel(1, 1, gallery.getSize(), 1));
//			jsColumns.setModel(new SpinnerNumberModel(1, 1, gallery.getSize(), 1));
//			jsGoToImage.setModel(new SpinnerNumberModel(1, 1, gallery.getSize(), 1));
//		}

		int TEXTWIDTH = 4;
		((JSpinner.NumberEditor) jsZoom.getEditor()).getTextField().setColumns(TEXTWIDTH);
		((JSpinner.NumberEditor) jsGoToImage.getEditor()).getTextField().setColumns(TEXTWIDTH);
		((JSpinner.NumberEditor) jsRows.getEditor()).getTextField().setColumns(TEXTWIDTH);
		((JSpinner.NumberEditor) jsColumns.getEditor()).getTextField().setColumns(TEXTWIDTH);
		
		initResliceButtonMenu();
		toolBar.add(reslicebt);
        searchbt = new JButton(XmippResource.getIcon("binocular.png"));
        searchbt.setEnabled(data.isTableMode());
        searchbt.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                try {
                    findReplace();
                } catch (Exception ex) {
                    Logger.getLogger(GalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        });
        toolBar.add(searchbt);
        plotbt = new JButton(XmippResource.getIcon("plot.png"));
        plotbt.setEnabled(data.isTableMode());
        plotbt.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                try {
                    plotColumns();
                } catch (Exception ex) {
                    Logger.getLogger(GalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        });
        toolBar.add(plotbt);
        chimerabt = new JButton(XmippResource.getIcon("chimera.png"));
        chimerabt.setEnabled(data.isVolumeMode());
        chimerabt.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                openWithChimera();
            }
        });
        toolBar.add(chimerabt);
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
						selectBlock((String)item);
				}

				@Override
				public Object getSelectedItem()
				{
					return data.getSelectedBlock();
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
		jcbVolumes.setModel(new VolumesComboBoxModel());
		cbPanel.add(jcbVolumes);
	}
        
        public void selectBlock(String block)
        {
            data.selectBlock(block);
            jcbVolumes.invalidate();
            try
            {
                    data.loadMd();
                    reloadTableData();
                    autoAdjustColumns(data.isAutoAdjust());
                    enableToolBarActions();
            }
            catch (Exception e)
            {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
            }
        }
        
        protected void enableToolBarActions()
        {
        	boolean hasRender = data.allowGallery();
        	boolean isCol = data.isColumnFormat();
        	
			btnChangeView.setEnabled(hasRender && isCol && data.renderImages());
			jsZoom.setEnabled(hasRender);
			jlZoom.setEnabled(hasRender);
			
			jsGoToImage.setEnabled(isCol && gallery.getSize() > 0);
			jlGoToImage.setEnabled(isCol);
			
        }
        
        public class VolumesComboBoxModel implements ComboBoxModel
        {
        	private HashMap<Integer, String> indexToVolume;
        	private HashMap<String, Integer> volumeToIndex;
        	private int size;
        	
        	public VolumesComboBoxModel()
        	{
        		String volume;
        		size = data.allowsVolumeMode()? data.size(): 0;
        		indexToVolume = new HashMap(size);
        		volumeToIndex = new HashMap(size);
        		for(int i = 0; i < size; i ++)
        		{
        			volume = data.getVolumeAt(i);
        			volume = getVolumeName(volume);
        			indexToVolume.put(i, volume);
        			volumeToIndex.put(volume, i);
        		}
        	}
            
        	

			@Override
			public int getSize()
			{
                    
                    return size;
			}

			@Override
			public Object getElementAt(int index)
			{
				return indexToVolume.get(index);
			}

			@Override
			public void setSelectedItem(Object anItem)
			{
				try
				{
					int index = volumeToIndex.get(anItem);
            		data.selectVolumeAt(index);
            		reloadTableData();
				}
                catch(Exception e)
                {
                	XmippDialog.showError(GalleryJFrame.this, e.getMessage());
                }
				
			}
                        
            public String getVolumeName(String volume)
            {

                    if(Filename.hasPrefix(volume))//is stack
                        return volume;
                    String base = Filename.getBaseName(volume);
                    int count = 0;
                    for(int i = 0; i < getSize(); i ++)
                        if(base.equals(Filename.getBaseName(data.getVolumeAt(i))))
                        {
                            count ++;
                            if (count > 1)// if volumes use same name you need to keep preffix
                            	break;
                        }
                    if(count == 1)
                        return base;
                    
                    return volume;//needs full path to differentiate volumes
            } 

			@Override
			public Object getSelectedItem()
			{
				return getVolumeName(data.getSelVolumeFile());
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
        }

	protected void updateVisibleCombos()
	{

		boolean showBlocks = data.getNumberOfBlocks() > 0;
		boolean showVols = data.size() > 1 && data.isVolumeMode();
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
			Integer rows = (Integer) jsRows.getValue();
			if (rows <= 0)
				rows = 1;
			else if (rows >= gallery.getSize())
				rows = gallery.getSize();
			if (jsRows.getValue() != rows)
				jsRows.setValue(rows);
			data.setModelDim(rows, null);
			gallery.setRows(rows);
		}
	}

	private void jsColumnsStateChanged(javax.swing.event.ChangeEvent evt)
	{// GEN-FIRST:event_jsColumnsStateChanged
		if (!isUpdating)
		{
			Integer columns = (Integer) jsColumns.getValue();
			if (columns <= 0)
				columns = 1;
			else if (columns >= gallery.getSize())
				columns = gallery.getSize();
			if (jsColumns.getValue() != columns)
				jsColumns.setValue(columns);
			data.setModelDim(null, columns);
			gallery.setColumns(columns);

		}
	}

	

	private void jsGoToImageStateChanged(javax.swing.event.ChangeEvent evt)
	{
		if (!isUpdating)
		{
			Integer oValue = ((Integer) jsGoToImage.getValue());
			Integer value = oValue;
			if (oValue <= 0)
			{
				value = 1;
				isUpdating = true;
			}
			else if (oValue >= gallery.getSize())
			{
				value = gallery.getSize();
				isUpdating = true;
			}
			if (oValue != value)
				jsGoToImage.setValue(value);
			if(oValue != 0 && oValue != gallery.getSelTo() + 1)//not me and not already selected
				goToImage(value - 1);
		}
		else
			isUpdating = false;
	}

	private void formComponentResized(java.awt.event.ComponentEvent evt)
	{
		width = getSize().width;
		if (!isUpdating && data.isAutoAdjust())
		{
			adjustColumns();
		}
	}

	public void selectItem(int row, int col)
	{
		if (row < 0 || row > table.getRowCount() - 1 || col < 0 || col > table.getColumnCount())
			return;
		
		if (gallery.data.isGalleryMode() && row * table.getColumnCount() + col + 1 > gallery.getSize())
		{
			Point coords = gallery.getCoords(gallery.getSize() - 1);
			row = coords.y;
			col = coords.x;
		}
		gallery.clearSelection();
		gallery.touchItem(row, col);
		makeVisible(row, col);

	}

	protected void tableMouseClicked(MouseEvent evt)
	{
		final Point p = evt.getPoint();
		int row = table.rowAtPoint(p);
		int col = table.columnAtPoint(p);
		col = table.convertColumnIndexToModel(col);
        int index = gallery.getIndex(row, col);
        if (!gallery.isValidIndex(index))
        	return;
        if (evt.getButton() == MouseEvent.BUTTON1 && evt.getClickCount() > 1)
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
                isUpdating = true;
                jsGoToImage.setValue(index + 1);
                isUpdating = false;
                // Ctrl adds items to selection, otherwise previous ones are removed.
                if (!evt.isControlDown() && !evt.isShiftDown())
                {
                        if(evt.getButton() == MouseEvent.BUTTON1)
                        {
                            gallery.clearSelection();
                            gallery.touchItem(row, col);
                            
                        }
                        else
                        {
                            if(!gallery.isSelected(row, col))
                                gallery.clearSelection();
                            gallery.touchItem(row, col, true);
                        }
                       
                }
                else
                {
                        if (evt.isShiftDown())
                                gallery.selectRange(previousSelectedRow, previousSelectedCol, row, col);
                        else if (evt.isControlDown())
                        {
                                gallery.touchItem(row, col);
                                jsGoToImage.setValue(gallery.getSelTo() + 1);
                        }
                }

                if (!evt.isShiftDown())
                {
                        previousSelectedRow = row;
                        previousSelectedCol = col;
                }

            final MouseEvent me = evt;
            if (evt.getButton() == MouseEvent.BUTTON3 && gallery.handleRightClick(row, col, jpopUpMenuTable))
            {
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

	private void autoAdjustColumns(boolean autoadjust)
	{
		setAutoAdjustColumns(autoadjust);
        if(autoadjust)
            data.setModelDim(null, null);
        else
        {
        	int rows = (Integer)jsRows.getValue();
        	int cols = (Integer)jsColumns.getValue();
        	if(rows > 0 && cols > 0)
        		data.setModelDim(rows, cols);
        }
		adjustColumns();
	}

	private void setResliceView(int view)
	{
		data.setResliceView(view);
		reloadTableData();
	}
        
    protected void plotColumns() {
        PlotJDialog dlg = new PlotJDialog(GalleryJFrame.this);
        dlg.showDialog();
    }

	protected class GalleryMenu extends XmippMenuBarCreator
	{

		protected QuickHelpJDialog quickhelpdlg;

		@Override
		protected void createItems() throws Exception
		{
			// File
			addItem(FILE, "File");
			addItem(FILE_OPEN, "Open ...", null, "control released O");
			addItem(FILE_OPENWITH_IJ, "Open with ImageJ", "ij.gif", "control released J");
			addItem(FILE_OPENWITH_CHIMERA, "Open with Chimera", "chimera.png", "control released H");
			addItem(FILE_OPENMICROGRAPHS, "Open colored particles");
			addItem(FILE_INFO, "File info ...");
			addExtraMenuItems();

			addSeparator(FILE);
			addItem(FILE_SAVE, "Save", "save.gif", "control released S");
			addItem(FILE_SAVEAS, "Save as", "save_as.gif");

            addItem(FILE_EXPORTIMAGES, "Export images ...", "export_wiz.gif");
			addItem(FILE_REFRESH, "Refresh", "refresh.gif", "released F5");
			addSeparator(FILE);
			addItem(FILE_EXIT, "Exit", null, "control released Q");
			// Display
			
                        
			addItem(DISPLAY, "Display");
            addItem(DISPLAY_APPLYGEO, "Apply geometry", null, "control released G");
			addItem(DISPLAY_WRAP, "Wrap", null, "control released W");
            addItem(DISPLAY_NORMALIZE, "Normalize", null, "control released N");
            addItem(DISPLAY_INVERTY, "Render positive Y up");
			addSeparator(DISPLAY);
			addItem(DISPLAY_RENDERIMAGES, "Render images", null, "control released R");
			
			addRenderImageColumnItems();
			addDisplayLabelItems();
			
			addItem(DISPLAY_RESLICE, "Reslice");
			for (int i = 0; i < ImageGeneric.VIEWS.length; ++i)
				addItem(DISPLAY_RESLICE_VIEWS[i], reslices[i]);
            addItem(DISPLAY_COLUMNS, "Columns ...", "columns.gif");
            
                                    
            addItem(TOOLS, "Tools");
			addItem(TOOLS_AVGSTD, "Avg & Std images");
			addItem(TOOLS_PCA, "PCA");
			addItem(TOOLS_FSC, "FSC");
			addItem(TOOLS_PLOT, "Plot", "plot.png");
			addItem(MD_FIND_REPLACE, "Find & Replace", "binocular.png", "control released F");
                        
			// Metadata operations
			addItem(METADATA, "Metadata");
                        
			
			addItem(MD_CLASSES, "Classes");
			addItem(MD_EDIT_COLS, "Edit labels", "edit.gif");
			addItem(MD_ADD_OBJECT, "Add new object", "new_object.gif");
			addItem(MD_REMOVE_DISABLED, "Remove disabled", "delete.gif");
			addItem(MD_REMOVE_SELECTION, "Remove selection");
			addItem(MD_SAVE_SELECTION, "Save state", "save.gif");
			addSeparator(METADATA);
			// Help
			addItem(HELP, "Help");
			addItem(HELP_ONLINE, "Online help", "online_help.gif");
			addItem(KEY_ASSIST, "Tips...", "bulb.png");
		}// function createItems
		
		//To insert extra items inside file menu
		public void addExtraMenuItems()
		{}

		public void update()
		{
                        
            boolean isscipion = data.isScipionInstance();
			boolean galMode = data.isGalleryMode();
			boolean volMode = !data.getSelVolumeFile().isEmpty();
			setItemEnabled(FILE_OPENWITH_CHIMERA, volMode || data.containsGeometryInfo("3D")|| data.containsGeometryInfo("Projection"));
			setItemEnabled(FILE_OPENMICROGRAPHS, data.hasMicrographParticles());
            setItemEnabled(FILE_EXPORTIMAGES, data.hasRenderLabel() && !volMode);
			setItemEnabled(FILE_SAVE, !volMode && !isscipion);
			setItemEnabled(FILE_SAVEAS, !volMode && !isscipion);
			
			setItemEnabled(DISPLAY_APPLYGEO, data.containsGeometryInfo());
			setItemEnabled(DISPLAY_WRAP, data.containsGeometryInfo() && data.useGeo());
			setItemSelected(DISPLAY_WRAP, data.containsGeometryInfo() && data.isWrap());
			setItemSelected(DISPLAY_APPLYGEO, data.useGeo());
            setItemSelected(DISPLAY_INVERTY, data.isInvertY());
			setItemEnabled(DISPLAY_RENDERIMAGES, data.isTableMode() && data.isColumnFormat());
			setItemSelected(DISPLAY_RENDERIMAGES, data.renderImages());
            setItemEnabled(DISPLAY_SHOWLABELS, gallery.showLabels());
			setItemEnabled(DISPLAY_RENDERIMAGE, galMode);
			for (int i = 0; i < ImageGeneric.VIEWS.length; ++i)
				setItemSelected(DISPLAY_RESLICE_VIEWS[i], (data.getResliceView() == ImageGeneric.VIEWS[i]));
			setItemEnabled(DISPLAY_COLUMNS, !galMode);
			setItemEnabled(DISPLAY_RESLICE, data.isVolumeMode());
            setItemSelected(DISPLAY_NORMALIZE, data.getNormalized());
			setItemEnabled(MD_CLASSES, data.isClassificationMd());
			setItemEnabled(TOOLS_PLOT, data.isTableMode());
			boolean isCol = data.isColumnFormat();
			boolean doStats = isCol && !volMode;
			setItemEnabled(TOOLS_AVGSTD, doStats);
			setItemEnabled(TOOLS_FSC, doStats);
			setItemEnabled(TOOLS_PCA, doStats);
			setItemEnabled(MD_ADD_OBJECT, isCol);
			setItemEnabled(MD_REMOVE_DISABLED, isCol);
			setItemEnabled(MD_REMOVE_SELECTION, isCol);
			setItemEnabled(MD_SAVE_SELECTION, isCol);
			setItemEnabled(MD_FIND_REPLACE, isCol && !galMode);
			reslicebt.setEnabled(data.isVolumeMode());
			chimerabt.setEnabled(volMode);
            setItemVisible(METADATA, !isscipion);
            addDisplayLabelItems();
            addRenderImageColumnItems();
		}// function update

		@Override
		protected void handleActionPerformed(ActionEvent evt)
		{
			String cmd = evt.getActionCommand();
			try
			{
                if (cmd.equals(DISPLAY_INVERTY))
				{
					gallery.setInvertY(getItemSelected(DISPLAY_INVERTY));
				}
                                else if (cmd.equals(DISPLAY_NORMALIZE))
				{
					gallery.setNormalized(getItemSelected(DISPLAY_NORMALIZE));
				}
				else if (cmd.equals(DISPLAY_APPLYGEO) || cmd.equals(DISPLAY_WRAP))
				{
					if (data.containsGeometryInfo())
					{
						((MetadataGalleryTableModel) gallery).setUseGeometry(getItemSelected(DISPLAY_APPLYGEO), getItemSelected(DISPLAY_WRAP));
						setItemEnabled(DISPLAY_WRAP, data.containsGeometryInfo() && data.useGeo());

					}
				}
				
				else if (cmd.equals(DISPLAY_RENDERIMAGES))
				{
                    if(gallery instanceof MetadataTableModel)
                    {
                        ((MetadataTableModel) gallery).setRenderImages(getItemSelected(DISPLAY_RENDERIMAGES));
                        setItemEnabled(DISPLAY_SHOWLABELS, gallery.showLabels());
                        btnChangeView.setEnabled(data.hasRenderLabel() && data.renderImages());
                    }
					makeVisible(gallery.getFirstSelectedIndex(), 0);
				}
				
				else if (cmd.equals(DISPLAY_COLUMNS))
				{
					ColumnsJDialog dialog = new ColumnsJDialog(GalleryJFrame.this);
					boolean result = dialog.showDialog();
					if (result)
					{
						List<ColumnInfo> columns = dialog.getColumnsResult();
						isUpdating = true;
						((MetadataGalleryTableModel) gallery).updateColumnInfo(columns);
						gallery.fireTableDataChanged();
						//setItemEnabled(DISPLAY_RENDERIMAGES, data.renderImages());
						// menu.enableRenderImages(data.globalRender);
						isUpdating = false;
					}
				}
				else if (cmd.equals(TOOLS_AVGSTD))
					runInBackground(Worker.STATS);
				else if (cmd.equals(TOOLS_PCA))
					runInBackground(Worker.PCA);
				else if (cmd.equals(TOOLS_FSC))
					runInBackground(Worker.FSC);
				else if (cmd.equals(FILE_OPEN))
				{
                    fc.setDialogTitle("Open");
                    fc.setApproveButtonToolTipText("File to open");
                    fc.setApproveButtonText("Open");
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
                else if (cmd.equals(FILE_EXPORTIMAGES))
				{
                	exportImages();
				}
				else if (cmd.equals(FILE_EXIT))
				{
					close();
				}
				else if (cmd.equals(FILE_OPENWITH_CHIMERA))
				{
                    openWithChimera();
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
                        XmippUtil.showImageJ(Tool.VIEWER);
						ImagePlusLoader loader = gallery.getImageLoader();
						ImagesWindowFactory.openXmippImageWindow(loader, data.parameters);
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
				else if (cmd.equals(TOOLS_PLOT))
				{
					plotColumns();
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
                                        if(gallery.hasSelection())
                                            data.saveSelection(gallery.getSelection());
				}
				else if (cmd.equals(MD_FIND_REPLACE))
				{
					findReplace();
				}
				else if (cmd.equals(MD_ADD_OBJECT))
				{
                                    boolean added = data.addObject();
                                    if(added)
                                        reloadMd();
				}
				else if (cmd.equals(HELP_ONLINE))
				{
					XmippWindowUtil.openURI("http://github.com/I2PC/scipion/wiki/ShowJ");
				}
				else if (cmd.equals(KEY_ASSIST))
				{
					if (quickhelpdlg == null)
						quickhelpdlg = new QuickHelpJDialog(GalleryJFrame.this, false, "Tips", getKeyAssist());
					quickhelpdlg.setVisible(true);

				}

			}
			catch (Exception e)
			{
				showException(e);
			}
		}// function handleActionPerformed
		
		
		protected void addRenderImageColumnItems()
		{                  
            JMenuItem mi = getItem(DISPLAY_RENDERIMAGE);
            if(mi == null)
                mi = addItem(DISPLAY_RENDERIMAGE, "Render label");
            else
                mi.removeAll();
			boolean rendercolumn = false;
			
            // Create the popup menu.
            String id, column;
            for(ColumnInfo ci: data.getColumns())
            {
                if(ci.render)
                {
                    column = ci.labelName;
                    id = String.format("Display.RenderImage.%s_rb", column.replace(".", ""));
                    mi = addItem(id, column);
                    mi.addActionListener(new RenderColumnActionListener());
                    if(data.getRenderColumn().toString().equals(column))
                            setItemSelected(id, true);
                    rendercolumn = true;
                }
            }
            setItemEnabled(DISPLAY_RENDERIMAGE, rendercolumn);

		}
                
        protected void addDisplayLabelItems()
		{                  
            JMenuItem mi = getItem(DISPLAY_SHOWLABELS);
			if(mi == null)
                mi = addItem(DISPLAY_SHOWLABELS, "Display label");
            else
                mi.removeAll();

            // Create the popup menu.
            String id, column;
            List<String> displayLabels = null;
            if(data.parameters.getDisplayLabels() != null)
                displayLabels = Arrays.asList(data.parameters.getDisplayLabels());
            for(ColumnInfo ci: data.getColumns())
            {
                if(ci.visible)
                {
                    column = ci.labelName;
                    id = String.format("Display.ShowLabel.%s_cb", column.replace(".", ""));
                    mi = addItem(id, column);
                    mi.addActionListener(new DisplayLabelActionListener());
                    
                    if(displayLabels!= null && displayLabels.contains(column))
                            setItemSelected(id, true);
                }   
            }

		}

        
                
        class DisplayLabelActionListener implements ActionListener
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				JCheckBoxMenuItem mi = (JCheckBoxMenuItem)e.getSource();
				String key = mi.getText();
                                
				data.setDisplayLabel(key, mi.isSelected());
				gallery.setShowLabels();
			}

		}
		
		
		class RenderColumnActionListener implements ActionListener
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				JRadioButtonMenuItem mi = (JRadioButtonMenuItem)e.getSource();
				String key = mi.getText();
				data.setRenderColumn(key);
                                
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
            setItemSelected(CTF_RECALCULATE, data.isRecalculateCTF(gallery.getIndex(row, col)));
			addSeparator();
			addItem(OPEN_IMAGES, "Open images");
			addItem(SAVE_IMAGES, "Save images", "save.gif");
			addItem(SET_CLASS, "Set class");
			addItem(SELECT, "Select");
			addItem(SELECT_ALL, "All", null, "control released A");
			addItem(SELECT_TOHERE, "To here");
			addItem(SELECT_FROMHERE, "From here");
                        addItem(INVERT_SELECT, "Invert selection");
                        if(data.parameters.objectCommands != null)
                            for(String cmd: data.parameters.objectCommands)
                                addItem(cmd + "_mi", cmd);
			initItems();
		}// function createItems

		public void show(Component cmpnt, Point location)
		{
            row = table.rowAtPoint(location);
			col = table.columnAtPoint(location);
			col = table.convertColumnIndexToModel(col);
            boolean isscipion = data.isScipionInstance();
			setItemVisible(SET_CLASS, data.isClassificationMd() && !isscipion);
			// This item visibility depends on current selection
			setItemVisible(SAVE_IMAGES, data.isClassificationMd() && gallery.getSelectionCount() > 0 && !isscipion);
			setItemVisible(OPEN_IMAGES, data.hasClasses() && gallery.getSelectionCount() == 1);
            setItemSelected(CTF_RECALCULATE, data.isRecalculateCTF(gallery.getIndex(row, col)));
			// Update menu items status depending on item.
			getPopupMenu().show(cmpnt, location.x, location.y);

		}// function show

		private void selectRange(int first, int last)
		{
			gallery.selectRange(first, last, true);
			jsGoToImage.setValue(gallery.getSelTo() + 1);
		}

	

		/** Set values to defaults */
		@Override
		public void initItems()
		{
			setItemVisible(OPEN, false);
			setItemVisible(OPEN_ASTEXT, false);
			setItemVisible(CTF_PROFILE, data.isCTFMd());
			setItemVisible(CTF_RECALCULATE, data.isCTFMd());
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
            else if (cmd.equals(INVERT_SELECT))
			{
				for (int i = 0; i < data.size(); i++)
                    gallery.setSelected(i, !gallery.isSelected(i));
                gallery.fireTableDataChanged();
                jsGoToImage.setValue(gallery.getSelTo() + 1);
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
				
                ColumnInfo ci = gallery.getColumn(row, col);
                if (ci.allowRender)
                    gallery.handleDoubleClick(row, col);
                else
                {
                    int index = gallery.getIndex(row, col);
                    String file = data.getValueFromCol(index, ci);
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
                                int index = gallery.getIndex(row, col);
				data.showCTF(true, index, gallery.getSelection(), ctfTasks);
			}
			else if (cmd.equals(CTF_RECALCULATE))
			{
                boolean isrecalculate = getItemSelected(CTF_RECALCULATE);
                int index = gallery.getIndex(row, col);
                if(isrecalculate && !data.isEnabled(index))
                    XmippDialog.showInfo(GalleryJFrame.this, "You must enable micrograph to recalculate its CTF");
                else
                {
                    if(isrecalculate)
                        data.showCTF(false, index, gallery.getSelection(), ctfTasks);
                    else
                        data.removeCTF(row);
                }
                                
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
                            if(gallery.hasSelection())
                            {
				File f = new File(data.getFileName());
				SaveImagesJDialog dialog = new SaveImagesJDialog(GalleryJFrame.this, f.getParent() + "/images_selection.xmd");
				dialog.showDialog();					
                            }
			}
                        
            else 
            {
                String objectCommand = cmd.replace("_mi", "");
                if (data.isObjectCmd(objectCommand))
                {
                    int index = gallery.getIndex(row, col);
                    data.runObjectCommand(index, objectCommand);
                }
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

		gallery.setRowBusy(row);
	}

	@Override
	public void setRowIdle(int row)
	{
		gallery.setRowIdle(row);

	}

	

	@Override
	public void done()
	{
		XmippDialog.showInfo(this, String.format("Calculating ctf: DONE"));
	}

	protected void saveMd() throws Exception
	{
		saveMd(dlgSave.getMdFilename(), false, dlgSave.isOverwrite(), true );
	}
        
      

	protected void saveMd(String path, boolean saveall, boolean isoverwrite, boolean reload) throws Exception
	{
		try
		{
			data.saveMd(path, saveall, isoverwrite);
            String file;
            if (path.contains("@"))
                    file = path.substring(path.lastIndexOf("@") + 1, path.length());
            else
            {
                    file = path;
                    path = getBlock() + "@" + file;
            }
            if(reload)
                {
                gallery.data.setFileName(file);
                if (path.contains("@"))
                        gallery.data.selectBlock(path.substring(0, path.lastIndexOf("@")));
                reloadFile(file, false);
                setGalleryTitle();
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
                data.saveAll(dlgSave.getMdFilename(), dlgSave.isOverwrite());
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

			
		}

	}
        
    private void exportImages() throws Exception
	{
		SwingUtilities.invokeLater(new Runnable() {

                    @Override
                    public void run() {
                        ExportImagesJDialog d = new ExportImagesJDialog(GalleryJFrame.this);
                    }
                });
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
		for (int i = 0; i < data.size(); i++)
			if (gallery.isSelected(i))
				extractframe.refreshActive(data.getId(i), data.isEnabled(i));

	}

	public void refreshActive(long id, boolean enabled)
	{
		for (int i = 0; i < data.size(); i++)
			if (id == data.getId(i))
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
		createModel(gallery.getSelection());
		reloadMd(changed);
		reloadCombos();
	}
	
	private void reloadCombos(){
		createCombos();
		updateVisibleCombos();
		jcbBlocks.setSelectedItem(gallery.data.getSelectedBlock());
		
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

	public Map<Object, Object> getKeyAssist()
	{
		Map<Object, Object> map = Collections.synchronizedMap(new LinkedHashMap<Object, Object>());
		map.put("Shift + scroll up/ctrl + P", "Zoom in if images displayed");
		map.put("Shift + scroll down/ctrl + M", "Zoom out if images displayed");
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
                map.put("Space Bar", "Enable/Disable selection");

		return map;
	}
        
        public TasksEngine getTasksEngine()
        {
            return ctfTasks;
        }
        
        
        
        protected void openChimera(String file, boolean link)
        {
            try
            {

                    String scipionHome = System.getenv().get("SCIPION_HOME");
                    if(scipionHome == null)
                    {
                        XmippDialog.showError(GalleryJFrame.this, "Scipion is not available");
                        return;
                    }
                    
                        
                    String run = String.format("python %s%2$spyworkflow%2$sapps%2$spw_chimera_client.py projector --input %3$s ", scipionHome, File.separator, file);
                    
                    if(data.parameters.getSamplingRate() != null)
                        run += String.format(" --samplingRate %s", data.parameters.getSamplingRate());
                    if(link)
                    {
                        int port = XmippUtil.findFreePort();
                        data.parameters.setChimeraPort(port);
                        run += "--showjPort " + port;
                    }
                    String output = XmippWindowUtil.executeCommand(run, false);
                    //System.out.println(output);

            }
            catch (Exception ex)
            {
                    ex.printStackTrace();
            }
        }
        
        protected void openWithChimera()
        {
        	if(data.containsGeometryInfo("3D") || data.containsGeometryInfo("Projection") )
            {
                fc.setApproveButtonText("Open");
                fc.setDialogTitle("Open with Chimera");
                fc.setApproveButtonToolTipText("Choose a chimera compatible file.");
                int result = fc.showOpenDialog(GalleryJFrame.this);
                if(result != XmippFileChooser.CANCEL_OPTION)
                {
                    String path = fc.getSelectedPath();
                    if(!Filename.isVolumeExt(path))
                    {
                        XmippDialog.showError(GalleryJFrame.this, XmippMessage.getFileTypeNotSupportedMsg(path));
                        return;
                    }
                    openChimera(path, true);
                }
            }
            else
                openChimera(data.getSelVolumeFile(), false);
        }
        
       
       
}// class JFrameGallery
