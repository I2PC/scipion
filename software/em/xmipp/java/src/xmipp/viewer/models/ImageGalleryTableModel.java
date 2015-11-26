/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

package xmipp.viewer.models;

import ij.ImagePlus;

import java.awt.Dimension;
import java.awt.Image;
import java.awt.Insets;
import java.awt.Point;
import java.util.ArrayList;

import javax.swing.JTable;
import javax.swing.border.Border;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.JTableHeader;

import xmipp.ij.commons.ImagePlusLoader;
import xmipp.jni.MDLabel;
import xmipp.jni.MDRow;
import xmipp.jni.MetaData;
import xmipp.utils.Cache;
import xmipp.utils.DEBUG;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.viewer.ImageDimension;
import xmipp.viewer.ImageItemRenderer;

public abstract class ImageGalleryTableModel extends AbstractTableModel {
    
    protected ArrayList<Integer> busyRows = new ArrayList<Integer>();

	private static final long serialVersionUID = 1L;
	// Store the number of rows and columns
	protected int rows, cols;
	// Total number of elements and last width
	protected int n, last_width;
	// Image original dimensions
	protected int image_width, image_height;
	// Thumbnails dimensions
	protected int thumb_width, thumb_height, font_height;
	// Table cell dimensions, this is redundant, but to avoid to be recalculated
	public Dimension cellDim = new Dimension();

	// Relation between real dimensions and thumbnails dimensions
	// scale = image_width / thumb_width
	protected float scale = (float) 1.;
	// protected int zoom = 100;
	// Cache class to reuse of already loaded items
	protected Cache<String, ImageItem> cache;
	// Filename
	// Hold gallery dimensions
	protected ImageDimension dimension;
	// Renderer to display images
	protected ImageItemRenderer renderer = new ImageItemRenderer();
	// Column model
	protected GalleryColumnModel columnModel;
	// Where to show labels
	// protected boolean showLabel = false;
	// Whether to autoadjust columns
	// Flags and variables to control global normalization
	protected boolean normalize_calculated = false;
	protected double normalize_min = Double.POSITIVE_INFINITY,
			normalize_max = Double.NEGATIVE_INFINITY;

	public GalleryData data; // information about the gallery

	public boolean adjustWidth = true; 
    protected boolean[] selection;
	
	// Initiazation function
	public ImageGalleryTableModel(GalleryData data, boolean[] selection) throws Exception {
		this.data = data;
		cols = 0;
		dimension = loadDimension();
		columnModel = createColumnModel();
		
		int zoom = 100;
		if(data.getZoom() != null)
			zoom = data.getZoom();
		setZoomValue(zoom);
		// This should be changed later after a call to
		// setColumns or adjustColumns
		if (cols < 1) {
			cols = 1;
			rows = n;
		}
		// DEBUG.printMessage(String.format("col: %d, rows: %d", cols, rows));
		//resizeCache(); NOw this is done when setZoomValue
		//might be that the number of items is the same and we keep selection, although it is a different metadata
        if(selection == null)
        	this.selection = new boolean[data.ids.length];
        else
        	this.selection = selection;
        
	}
	
	
	
	public int getImageWidth()
	{
		return image_width;
	}
	
	public int getImageHeight()
	{
		return image_height;
	}
	
	// Load initial dimensions
	protected abstract ImageDimension loadDimension() throws Exception;

	/**
	 * Function to create Cache of items calculate the dimension of the cache
	 * depending on images size and available memory
	 */
	protected void resizeCache() {
		int limit = Cache.getLimit(thumb_width, thumb_height);
		if(cache == null)
			cache = new Cache(limit);
		else
			cache.resize(limit > 0 ? limit : 1);
	}

	@Override
	public int getColumnCount() {
		return cols;
	}

	@Override
	public int getRowCount() {
		return rows;
	}

	@Override
	public String getColumnName(int column) {
		return String.format("%d", column + 1);
	}

	@Override
	public Class getColumnClass(int c) {
		return ImageItem.class;
	}

	@Override
	public Object getValueAt(int row, int col) {
		int index = getIndex(row, col);
		if (isValidIndex(index)) {
			try {
				String key = getItemKey(index);
				ImageItem item;
				// If the element is on cache, just return it
				if (cache.containsKey(key))
					item = cache.get(key);
				else {
//					// If not, create the item and store it for future
					item = createItem(index, key, getColumn(row, col));
					cache.put(key, item);
				}
				setupItem(item);
				return item;
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		return null;
	}

	/** Function to force the refresh of some item */
	public void refreshAt(int row, int col) {
		
		int index = getIndex(row, col);

		if (isValidIndex(index)) {
			try {
				String key = getItemKey(index, data.getLabelFromCol(col));
				// Clear cache entry
				cache.remove(key);
				fireTableRowsUpdated(row, row);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}



	protected void setupItem(ImageItem item) {
		ImagePlus imp = item.getImagePlus();
		if (imp != null && imp.getProcessor() != null) { // When image is missing this will be null
			if (data.normalize)
				imp.getProcessor().setMinAndMax(normalize_min, normalize_max);
			else
				imp.getProcessor().resetMinAndMax();
			imp.updateImage();
		}
	}

	public void setRows(int rows) {
		if (rows != this.rows) {
			this.rows = rows;
			cols = n / rows + (n % rows == 0 ? 0 : 1);
			fireTableStructureChanged();
		}
	}

	protected void setColumnsValue(int cols) {
		this.cols = cols;
		rows = n / cols + (n % cols == 0 ? 0 : 1);
		fireTableStructureChanged();
	}

	public void setColumns(int cols) {
		if (cols != this.cols) {
			setColumnsValue(cols);
		}
	}

	// Return True if columns were changed
	public boolean adjustColumn(int width) {
		last_width = width;
		int new_cols = Math.max(Math.min(width / cellDim.width, n), 1);
		if (new_cols != cols) {
			setColumnsValue(new_cols);
			return true;
		}
		return false;
	}

	public int getSize() {
		return n;
	}

	public Dimension getCellSize() {
		return cellDim;
	}// function getCellSize

	/**
	 * This function will be used to calculate all size need to render the image
	 * cell
	 * 
	 * @param width
	 * @param height
	 */
	protected void calculateCellSize() {
		thumb_width = (int) (image_width * scale);
		thumb_height = (int) (image_height * scale);

		font_height = 0;
		if (data.isDisplayLabel()) {
			font_height = renderer.getFontMetrics(renderer.getFont())
					.getHeight();
			font_height += renderer.getIconTextGap(); // Adds the extra gap.
			// font_height -= table.getIntercellSpacing().height; // Removes
		}

		int borderHeight = 0, borderWidth = 0;
		Border border = renderer.getBorder();

		if (border != null) {
			Insets insets = renderer.getBorder().getBorderInsets(renderer);
			borderWidth = insets.left + insets.right;
			borderHeight = insets.bottom + insets.top;
		}
		cellDim.setSize(thumb_width + 1 * borderWidth, thumb_height + 1
				* borderHeight + font_height);
		adjustColumnsWidth();
	}

	/**
	 * This method will be used to set renderers and other table configurations
	 */
	public void setupTable(JTable table) {
		table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		table.setDefaultRenderer(ImageItem.class, renderer);
	}

	/** Update the table selection according with data selection */
	public void updateTableSelection(JTable table) {
		// For now do nothing, overrided in MetadataTable
	}

	/** Ajust the width of columns and headers */
	protected void adjustColumnsWidth() {
		// renderer.setPreferredSize(cellDim);
		columnModel.setWidth(cellDim.width);
	}

	/** Return the cell renderer to be used to draw images */
	public ImageItemRenderer getRenderer() {
		return renderer;
	}

	/** Return the column model to be used with this table model */
	public GalleryColumnModel getColumnModel() {
		return columnModel;
	}
	
	public JTableHeader getTableHeaderModel(){
		return new JTableHeader(columnModel);
	}

	/** Internal method to create the column model */
	protected GalleryColumnModel createColumnModel() {
		return new GalleryColumnModel(cellDim.width);
	}

	protected void setZoomValue(int z) {
		data.setZoom(z);
		scale = (float) (z / 100.0);
		calculateCellSize();
		resizeCache();
	}

	public void setZoom(int z) {
		
		if (data.getZoom() != null || data.getZoom() != z) {
			setZoomValue(z);
			fireTableDataChanged();
			if (data.isAutoAdjust())
				adjustColumn(last_width);
		}
	}

	/**
	 * Calculate index base on x position and y position
	 * 
	 * @param x
	 *            x position on table
	 * @param y
	 *            y position on table
	 * @return index of the element
	 */
	public int getIndex(int row, int col) {
                
		return row * cols + col; 
	}
        
       
	
	public long getId(int row, int col) {
		int index = getIndex(row, col);
		return data.ids[index];
		
	}
	
	/** Check if the index is inside bounds */
	public boolean isValidIndex(int index){
		return index >= 0 && index < n;
	}

	/** Return the label of the item at this position */
	public abstract String getLabel(int row, int col);

	/**
	 * Return the row and col given an index
	 * 
	 * @param index
	 * @return
	 */
	public Point getCoords(int index) {
                Point p = new Point();
		p.y = index / cols;
		p.x = index % cols;
		return p;
	}

	/** Functions to handle selections */

	/** Set enable/disable to selected items */
	public void setSelectionEnabled(boolean value) {
		try {
                    int from = getSelFrom(), to = getSelTo();
			for (int i = from; i <= to; i++)
				if (selection[i])
					data.setEnabled(i, value);
                                
		    if(from != -1)
		    {
		        from = getCoords(from).y;
		        to = getCoords(to).y;
		        fireTableRowsUpdated(from, to);
		    }
	            
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	
	
	public int getFirstSelectedIndex(){
		return getSelFrom();
	}

	/** Select a range of elements given the indexes */
	public void selectRange(int first_index, int last_index, boolean value) {
		for (int i = first_index; i <= last_index; ++i)
			setSelected(i, value);
                int from = getCoords(first_index).y;
                int to = getCoords(last_index).y;
		fireTableRowsUpdated(from, to);
	}

	/** Select a range of elements given the coordinates */
	public void selectRange(int first_row, int first_col, int last_row, int last_col) {
		int i1 = getIndex(first_row, first_col);
		int i2 = getIndex(last_row, last_col);
		int i = Math.min(i1, i2);
		i2 = Math.max(i1, i2);
		for (; i <= i2; ++i)
			setSelected(i, true);
		fireTableRowsUpdated(first_row, last_row);
	}

	/** Set the class of a given selection of elements. */
	public void setSelectionClass(int classNumber) {
		ClassInfo cli = data.getClassInfo(classNumber);
		for (int i = getSelFrom(); i <= getSelTo(); ++i)
			if (selection[i]) {
				data.setItemClass(i, cli);
			}
		clearSelection();
		fireTableDataChanged();
	}

	/** Remove a class */
	public void removeClass(int classNumber) {
		data.removeClass(classNumber);
		fireTableDataChanged();
	}

	/** Set the selection state of an element give row and col */
	public void touchItem(int row, int col) {
        setSelected(row, col, !isSelected(row, col));
        adjustWidth = false;
        fireTableCellUpdated(row, col);
	}
        
        /** Set the selection state of an element give row and col */
	public void touchItem(int row, int col, boolean isselected) {
		
        setSelected(row, col, isselected);
        adjustWidth = false;
        fireTableCellUpdated(row, col);	
	}
	
	
	/**
	 * Goto and select specified item, if there is a selection it will be
	 * cleared
	 */
	public void gotoItem(int i) {
		clearSelection();
		setSelected(i, !isSelected(i));
		Point coords = getCoords(i);
		fireTableCellUpdated(coords.y, coords.x);
	}

	

	/** Normalization utils */
	public void setNormalized(boolean normalize) {

		if (normalize != data.normalize) {
			data.normalize = normalize;
			if (normalize)
				calculateMinAndMax();
			fireTableDataChanged();
		}
	}

	

	/** Calculate min and max if needed for normalization */
	public void calculateMinAndMax() {
		if (!normalize_calculated) {
			double[] mm = getMinAndMax();
			normalize_min = mm[0];
			normalize_max = mm[1];
			normalize_calculated = true;
		}
	}

	/**
	 * This function will be called when a right click is fired from the UI, the
	 * col and row of cell where the event was produced and also the JPopupMenu
	 * This function will serve to personalize the menu before showing or event
	 * not show at all if the return is false
	 */
	public abstract boolean handleRightClick(int row, int col,
		XmippPopupMenuCreator xpopup);

	/**
	 * This function will be called when a double click is performed under the
	 * element at specified row and column
	 */
	public abstract boolean handleDoubleClick(int row, int col) ;

	/** Whether to display the labels */
	public void setShowLabels() {
		
		calculateCellSize();
		fireTableStructureChanged();
		
	}

	

	/** Retrieve the mininum and maximum of data */
	protected abstract double[] getMinAndMax();

	// /** Retrieve the image filename, return null if error */
	// protected abstract String getImageFilenameAt(int row, int col);

	/**
	 * Function to create the key of the item knowing the item index
	 */
	public abstract String getItemKey(int index) throws Exception;
	
	/**
	 * Return a key string using label
	 */
	public String getItemKey(int index, int label) throws Exception {
		Object value = data.getValueFromLabel(index, label);
		
		String format = value + "_%d_%d(%d,%d)";
		if (data.useGeo)
			format += "_geo";
		if (data.wrap)
			format += "_wrap";
		// String key = String.format(format, thumb_width, thumb_height);
		// DEBUG.printMessage(String.format("key: %s", key));
        String key = String.format(format, index, label, thumb_width, thumb_height);
		return key;
	}

	
	/**
	 * Return the main title to be used on windows
	 */
	public abstract String getTitle();

	/**
	 * Return the ImagePlusLoader subclass that handle the reading
	 * of the gallery representation as an ImagePlus
	 */
	public abstract ImagePlusLoader getImageLoader();

	/**
	 * Function to create the image item. this should be implemented in
	 * subclasses of ImageGallery
	 * 
	 * @param index
	 *            The index of the item to create
	 * @return The item created
	 */
	protected abstract ImageItem createItem(int index, String key, ColumnInfo ci)
			throws Exception;

        public boolean isSelected(int row, int col) {
            int index = getIndex(row, col);
            return isSelected(index);
        }

        public abstract boolean showLabels() ;

        public void setInvertY(boolean itemSelected) {
            data.inverty = itemSelected;
            cache.clear();
            fireTableDataChanged();
        }

    

	/**
	 * This class will contains basic info to be used for image rendering. It
	 * will contains an ImagePlus, label and some other useful info.
	 */
	public class ImageItem {

		protected ImagePlus image = null;
		protected int index;

		/**
		 * First argument is the gallery to which this item belongs Constructor
		 * of ImageItem */
		public ImageItem(int row, int col) {
			index = ImageGalleryTableModel.this.getIndex(row, col);
		}

		public ImageItem(int index, ImagePlus imp) {
			this.index = index;
			this.image = imp;
		}

		public int getIndex() {
			return index;
		}

		public ImagePlus getImagePlus() {
			return image;
		}

		public void setImagePlus(ImagePlus value) {
			image = value;
		}

		public Image getImage() {
			if (image != null)
				return image.getImage();
			return null;
		}

		public boolean getShowLabel() {
			return data.isDisplayLabel();
		}

		public Dimension getCellDim() {
			return cellDim;
		}

		public String getLabel() {
			Point coords = getCoords(index);
			return ImageGalleryTableModel.this.getLabel(coords.y, coords.x);
		}

		public String getKey() {
			try {
				return getItemKey(index);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return null;
		}

		public boolean isSelected() {
			return ImageGalleryTableModel.this.isSelected(index);
		}

		public boolean isEnabled() {
			return data.isEnabled(index);
		}

		public boolean isBusy() {
			return ImageGalleryTableModel.this.isBusy(index, 0);
		}

		public ClassInfo getClassInfo() {
			return data.getItemClassInfo(index);
		}

		@Override
		public String toString() {
			return image != null ? image.getFileInfo().fileName : null;
		}
	}// class ImageItem

	public int getRows()
	{
		return rows;
	}


	public int getColumns()
	{
		return cols;
	}

        
                /** Check if the item is busy */
	public boolean isBusy(int row, int col){
		return busyRows.contains(row);
	}
    
	public void setRowBusy(int row) {
		DEBUG.printFormat("setting busy row: %d", row);
		busyRows.add(new Integer(row));
		fireTableRowsUpdated(row, row);
	}

	public void setRowIdle(int row) {
		DEBUG.printFormat("setting idle row: %d", row);
		busyRows.remove(new Integer(row));	
		long objId = data.ids[row];
		String psdFile = data.md.getValueString(MDLabel.MDL_PSD, objId);
		String sortFn = psdFile.replace(".psd", ".xmd");
		MetaData mdRow = new MetaData(sortFn);
		MDRow row2 = new MDRow();
		mdRow.getRow(row2, mdRow.firstObject());
		data.setRow(row2, objId);
		mdRow.destroy();
		refreshRow(row);
	}
	
	/** Function to force the refresh of some row */
	public void refreshRow(int row) {
		try {
			int n = getColumnCount();
			int index;
			String key;
			ColumnInfo ci;
			for (int col = 0; col < n; ++col) {
				ci = data.getColumnInfo(col);
				if (ci.allowRender){
					key = getItemKey(row, ci.label);
					DEBUG.printFormat("Removing item: %s from cache", key);
					cache.remove(key);
				}
			}
			fireTableRowsUpdated(row, row);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
        
        public void clearSelection() {
            for (int i = 0; i < selection.length; ++i) {
                selection[i] = false;
            }

        }
        
        /**
     * Return the number of selected elements
     */
    public int getSelectionCount() {
        int count = 0;

        if (!data.isVolumeMode() && hasSelection()) {
            for (int i = 0; i < selection.length; i++) {
                if (selection[i]) {
                    count++;
                }
            }
        }
        return count;
    }
        
    public boolean[] getSelection()
    {
        return selection;
    }
    
    public void setSelected(int index, boolean isselected) {
        
        
        selection[index] = isselected;
       
        if (data.isVolumeMd && data.isTableMode())
            data.selectedVolFn = isselected? data.getVolumeAt(index): data.getVolumeAt(0);

    }
    
    public boolean isSelected(int index) {
        
        return selection[index];
    }
    
    public int getSelFrom()
    {
    	for(int i = 0; i < selection.length; i ++)
    	{
    		if(selection[i])
	    		return i;
    	}
    	return -1;
    }

    public int getSelTo()
    {
    	int selto = -1;
    	for(int i = 0; i < selection.length; i ++)
    	{
    		if(selection[i])
    			selto = i;
    	}
        return selto;
    }
    
      public boolean hasSelection() {
        if(getSelFrom() == -1)
            return false;
        return true;
    }
      
    public void setSelected(int row, int col, boolean b) {
          setSelected(getIndex(row, col), b);
    }

	public ColumnInfo getColumn(int row, int col)
	{
		return data.ciFirstRender;
	}

}
