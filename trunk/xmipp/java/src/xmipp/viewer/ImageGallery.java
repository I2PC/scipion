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

package xmipp.viewer;

import ij.ImagePlus;

import java.awt.Dimension;
import java.awt.Insets;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.List;

import javax.swing.border.Border;
import javax.swing.table.AbstractTableModel;

import xmipp.utils.Cache;
import xmipp.utils.DEBUG;

public abstract class ImageGallery extends AbstractTableModel {

	private static final long serialVersionUID = 1L;
	// Store the number of rows and columns
	protected int rows, cols;
	// Total number of elements and last width
	protected int n, last_width;
	// Image original dimensions
	protected int image_width, image_height;
	// Thumbnails dimensions
	protected int thumb_width, thumb_height;
	// Table cell dimensions, this is redundant, but to avoid to be recalculated
	protected Dimension cellDim = new Dimension();

	// Relation between real dimensions and thumbnails dimensions
	// scale = image_width / thumb_width
	protected float scale = (float) 1.;
	protected int zoom = 100;
	// Cache class to reuse of already loaded items
	protected Cache<String, ImageItem> cache = new Cache<String, ImageItem>();
	// Filename
	protected String filename;
	// Hold gallery dimensions
	protected ImageDimension dimension;
	// Renderer to display images
	protected ImageItemRenderer renderer = new ImageItemRenderer();
	// Where to show labels
	protected boolean showLabel = false;
	// Whether to autoadjust columns
	protected boolean adjustColumns = false;
	// Store the selection state for each item
	protected boolean[] selection;
	// Flags and variables to control global normalization
	protected boolean normalize = false, normalize_calculated = false;
	protected double normalize_min = Double.POSITIVE_INFINITY,
			normalize_max = Double.NEGATIVE_INFINITY;

	// Initiazation function
	public ImageGallery(String fn, int zoom) throws Exception {
		filename = fn;
		dimension = loadDimension();
		// Zdim will always be used as number of elements to display
		n = dimension.getZDim();
		selection = new boolean[n];
		image_width = dimension.getXDim();
		image_height = dimension.getYDim();
		setZoomValue(zoom);
		// This should be changed later after a call to
		// setColumns or adjustColumns
		cols = 1;
		rows = n;
		// DEBUG.printMessage(String.format("col: %d, rows: %d", cols, rows));
		resizeCache();
	}
	
	// Load initial dimensions
	protected abstract ImageDimension loadDimension() throws Exception;

	/**
	 * Function to create Cache of items calculate the dimension of the cache
	 * depending on images size and available memory
	 */
	protected void resizeCache() {
		int imageSize = dimension.getXDim() * dimension.getYDim()
				* Cache.MAXPXSIZE;
		int elements = Cache.MEMORY_SIZE / imageSize;
		cache.resize(elements > 0 ? elements : 1);
	}

	@Override
	public int getColumnCount() {
		// DEBUG.printMessage(String.format("count requested: cols:%d", cols));
		return cols;
	}

	@Override
	public int getRowCount() {
		// DEBUG.printMessage(String.format("count requested: rows:%d", rows));
		// DEBUG.printStackTrace();
		return rows;
	}

	@Override
	public String getColumnName(int column) {
		return "" + (column + 1);
	}

	@Override
	public Class getColumnClass(int c) {
		return ImageItem.class;
	}

	@Override
	public Object getValueAt(int row, int col) {
		int index = getIndex(row, col);

		if (index >= n)
			return null;

		String key = getItemKey(index);
		// DEBUG.printMessage(String.format("with key %s", key));
		ImageItem item;
		// If the element is on cache, just return it
		if (cache.containsKey(key))
			item = cache.get(key);
		else {
			// If not, create the item and store it for future
			item = createItem(index, key);
			cache.put(key, item);
		}
		item.isSelected = selection[index];
		item.showLabel = showLabel;
		ImagePlus imp = item.getImage();
		if (normalize)
			imp.getProcessor().setMinAndMax(normalize_min, normalize_max);
		else
			imp.getProcessor().resetMinAndMax();
		imp.updateImage();
		
		return item;
	}

	public void setRows(int rows) {
		adjustColumns = false;
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
		// fireTableDataChanged();
	}

	public void setColumns(int cols) {
		adjustColumns = false;
		if (cols != this.cols) {
			setColumnsValue(cols);
		}
	}

	public void adjustColumn(int width) {
		last_width = width;
		adjustColumns = true;
		int new_cols = width / cellDim.width;
		if (new_cols != cols) {
			setColumnsValue(new_cols);
		}
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
	private void calculateCellSize() {
		thumb_width = (int) (image_width * scale);
		thumb_height = (int) (image_height * scale);

		int font_height = 0;
		if (showLabel) {
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
		renderer.setPreferredSize(cellDim);
	}

	public ImageItemRenderer getRenderer() {
		return renderer;
	}

	protected void setZoomValue(int z) {
		zoom = z;
		scale = (float) (zoom / 100.0);
		calculateCellSize();
	}

	public void setZoom(int z) {
		if (zoom != z) {
			setZoomValue(z);
			fireTableDataChanged();
			if (adjustColumns)
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

	/**
	 * Return the x and y given an index
	 * 
	 * @param index
	 * @return
	 */
	public int[] getCoords(int index) {
		int[] coords = new int[2];
		coords[0] = index % cols;
		coords[1] = index / cols;
		return coords;
	}

	/** Functions to handle selections */

	/** Clear selection list */
	public void clearSelection() {
		for (int i = 0; i < n; ++i)
			selection[i] = false;
	}

	/** Select a range of elements */
	public void touchRange(int first_row, int first_col, int last_row,
			int last_col) {
		int i1 = getIndex(first_row, first_col);
		int i2 = getIndex(last_row, last_col);
		int i = Math.min(i1, i2);
		i2 = Math.max(i1, i2);
		for (; i <= i2; ++i)
			selection[i] = !selection[i];
		;
		fireTableDataChanged();
	}

	/** Set the selection state of an element */
	public void touchItem(int row, int col) {
		int i = getIndex(row, col);
		selection[i] = !selection[i];
		fireTableCellUpdated(row, col);
	}

	/** Return the number of selected elements */
	public int getSelectedCount() {
		int count = 0;
		for (int i = 0; i < n; ++i)
			if (selection[i])
				++count;
		return count;
	}

	/** Normalization utils */
	public void setNormalized(boolean normalize) {

		if (normalize != this.normalize) {
			this.normalize = normalize;
			if (normalize)
				calculateMinAndMax();
			fireTableDataChanged();
		}
	}
	
	/** Calculate min and max if needed for normalization */
	public void calculateMinAndMax(){
		if (!normalize_calculated) {
			double[] mm = getMinAndMax();
			normalize_min = mm[0];
			normalize_max = mm[1];
			normalize_calculated = true;
		}
	}

	/** Whether to display the labels */
	public void setShowLabels(boolean value){
		if (showLabel != value){
			showLabel = value;
			calculateCellSize();
			fireTableDataChanged();
		}
	}
	
	/** Retrieve the mininum and maximum of data */
	protected abstract double[] getMinAndMax();
	

	/**
	 * Function to create the key of the item knowing the item index
	 */
	protected abstract String getItemKey(int index);

	/**
	 * Function to create the image item. this should be implemented in
	 * subclasses of ImageGallery
	 * 
	 * @param index
	 *            The index of the item to create
	 * @return The item created
	 */
	protected abstract ImageItem createItem(int index, String key);
}
