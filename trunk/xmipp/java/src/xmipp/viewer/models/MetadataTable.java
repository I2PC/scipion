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

import java.awt.Component;

import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.viewer.ImageItem;
import xmipp.viewer.TestRenderer;

public class MetadataTable extends MetadataGallery {
	private static final long serialVersionUID = 1L;

	public MetadataTable(GalleryData data) throws Exception {
		super(data);
		cols = visibleLabels.size();
	}

	@Override
	public String getColumnName(int column) {
		return visibleLabels.get(column).getLabelName();
	}

	@Override
	public Class getColumnClass(int column) {
		ColumnInfo ci = visibleLabels.get(column);
		if (ci.render)
			return ImageItem.class;
		else if (ci.getLabel() == MDLabel.MDL_ENABLED)
			return Boolean.class;
		try {
			return MetaData.getLabelClass(ci.getLabel());
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	@Override
	public int getIndex(int row, int col) {
		return row;
	}
	
	@Override
	public Object getValueAt(int row, int column) {
		// DEBUG.printMessage(String.format("MetadataTable.getValueAt(%d, %d)",
		// row, column));
		try {
			ColumnInfo ci = visibleLabels.get(column);
			if (ci.render) {
				String key = getItemKey(row, ci.getLabel());
				ImageItem item;
				// If the element is on cache, just return it
				if (cache.containsKey(key))
					item = cache.get(key);
				else {
					// If not, create the item and store it for future
					item = createImageItem(row, ci.getLabel(), ci.getLabel(), key);
					cache.put(key, item);
				}
				setupItem(item, row);
				return item;
			}
			int label = ci.getLabel();
			long id = data.ids[row];
			int type = MetaData.getLabelType(label);
			MetaData md = data.md;
			switch (type) {
			case MetaData.LABEL_INT:
				int value = md.getValueInt(label, id);
				// treat special case of MDL_ENABLED
				if (label == MDLabel.MDL_ENABLED)
					return (value > 0);
				return value;
			case MetaData.LABEL_BOOL:
				return md.getValueBoolean(label, id);
			case MetaData.LABEL_FLOAT:
			case MetaData.LABEL_DOUBLE:
				return md.getValueDouble(label, id);
			case MetaData.LABEL_LONG:
				return md.getValueLong(label, id);
			case MetaData.LABEL_STRING:
			case MetaData.LABEL_VECTOR:
			case MetaData.LABEL_VECTOR_LONG:
				return md.getValueString(ci.getLabel(), data.ids[row]);

			}
			return null;

		} catch (Exception e) {
			e.printStackTrace();
		}

		return null;
	}//function getValueAt

	@Override
	public void setValueAt(Object value, int row, int column) {
		try {
			ColumnInfo ci = visibleLabels.get(column);
			long id = data.ids[row];
			setMdValueAt(value, row, column, ci, id);
		} catch (Exception e) {
			e.printStackTrace();
		}

	}// function setValueAt
	
	/** Helper function to set the value of a cell 
	 * @throws Exception */
	protected void setMdValueAt(Object value, int row, int column,
			ColumnInfo ci, long id) throws Exception{
		if (!ci.render) {
			int label = ci.getLabel();
			int type = MetaData.getLabelType(label);
			MetaData md = data.md;
			switch (type) {
			case MetaData.LABEL_BOOL:
				md.setValueBoolean(label, (Boolean) value, id);
				break;
			case MetaData.LABEL_INT:
				if (label == MDLabel.MDL_ENABLED) {
					md.setEnabled((Boolean)value, id);
				} else
					md.setValueInt(label, ((Integer) value).intValue(), id);
				break;
			case MetaData.LABEL_FLOAT:
			case MetaData.LABEL_DOUBLE:
				md.setValueDouble(label, ((Double) value).doubleValue(), id);
				break;
			case MetaData.LABEL_LONG:
				md.setValueInt(label, ((Integer) value).intValue(), id);
				break;
			case MetaData.LABEL_STRING:
			case MetaData.LABEL_VECTOR:
			case MetaData.LABEL_VECTOR_LONG:
				// TODO: Implement a better editor for vectors
				md.setValueString(label, value.toString(), id);
				break;
			}
			fireTableRowsUpdated(row, row);
		}
	}

	@Override
	public boolean isCellEditable(int row, int column) {
		ColumnInfo ci = visibleLabels.get(column);
		return ci.allowEdit && !ci.render;
	}

	@Override
	public String getImageFilenameAt(int row, int col){
		ColumnInfo ci = visibleLabels.get(col);
		return (ci.render && data.isImageFile(ci))
				? data.getValueFromCol(row, ci) : null;
	}
	
	@Override
	public void setColumns(int cols) {

	}

	@Override
	public boolean adjustColumn(int width) {
		return false;
	}

	@Override
	/** Whether to display the labels */
	public void setRenderImages(boolean value) {
		boolean changed = false;
		for (ColumnInfo ci : data.labels)
			if (ci.allowRender && ci.render != value) {
				ci.render = value;
				changed = true;
			}
		if (changed) {
			data.globalRender = value;
			calculateCellSize();
			fireTableDataChanged();
		}
	}

	@Override
	protected void calculateCellSize() {
		//DEBUG.printMessage(String.format("MetadataTable:calculateSize"));
		if (data.globalRender) {
			super.calculateCellSize();
//			DEBUG.printMessage(String.format(
//					"MetadataTable:calculateSize w:%d, h:%d", cellDim.width,
//					cellDim.height));

		} else {
			int font_height;
			font_height = renderer.getFontMetrics(renderer.getFont())
					.getHeight();
			font_height += renderer.getIconTextGap(); // Adds the extra gap.
			cellDim.setSize(100, font_height + 5);
		}
	}

	@Override
	public void setupTable(JTable table) {
		//table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
		table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		table.setDefaultRenderer(ImageItem.class, renderer);
		table.setDefaultRenderer(Double.class, new TestRenderer());
		table.setAutoCreateRowSorter(true);
		updateTableSelection(table);
	}
	
	/** Update the table selection according with data selection */
	@Override
	public void updateTableSelection(JTable table){
		table.clearSelection();
		for (int i = 0; i < n; ++i)
			if (data.selection[i]) {
				table.addRowSelectionInterval(i, i);
			}
	}
	
	@Override
	public boolean handleRightClick(int row, int col, XmippPopupMenuCreator xpopup) {
		xpopup.initItems();
		if (data.isFile(col)){
			xpopup.setItemVisible(XmippPopupMenuCreator.OPEN, true);
			if (!data.isImageFile(col))
				xpopup.setItemVisible(XmippPopupMenuCreator.OPEN_ASTEXT, true);
		}
		return true;
	}

	@Override
	/** Return the column model to be used with this table model */
	public GalleryColumnModel createColumnModel() {
		return new MetadataColumnModel();
	}
	

	public class MetadataColumnModel extends GalleryColumnModel {
		public MetadataColumnModel() {
			super(0);
		}

		@Override
		public void adjustColumnsWidth(JTable table) {
			try {
				calculateCellSize();
				// String[] row = md.getRowValues(data.ids[0]);
				// table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
				int width = 0;
				TableCellRenderer rend;
				Component comp;
				for (int i = 0; i < visibleLabels.size(); ++i) {
					ColumnInfo col = visibleLabels.get(i);
					//Calculate width of the cell 
					if (col.render) {
						width = cellDim.width;
					} else {
						rend = table.getCellRenderer(0, i);
						comp = rend.getTableCellRendererComponent(
								table, getValueAt(0, i), false, false, 0, 0);
						width = comp.getPreferredSize().width + 10;
					}
					//Calculate width of the header
					TableColumn tc = getColumn(i);
					rend = tc.getHeaderRenderer();
					if (rend == null)
						rend = table.getTableHeader().getDefaultRenderer();
					Object value = tc.getHeaderValue();
					comp = rend.getTableCellRendererComponent(table, value, false, false, 0, i);
					//Take max width
					width = Math.max(width, comp.getPreferredSize().width);
					getColumn(i).setPreferredWidth(width);
//					DEBUG.printMessage(String.format("col: %d, width: %d", i,
//							width));
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}//function adjustColumnsWidth

	}//class MetadataColumnModel

}//class MetadataTable
