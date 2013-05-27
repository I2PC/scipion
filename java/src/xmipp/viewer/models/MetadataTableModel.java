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
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JTable;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;

import xmipp.ij.commons.ImagePlusLoader;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.viewer.FloatRenderer;
import xmipp.viewer.models.MetadataGalleryTableModel.MdRowImageLoader;
import xmipp.viewer.windows.ImagesWindowFactory;


public class MetadataTableModel extends MetadataGalleryTableModel
{
	private static final long serialVersionUID = 1L;

	int sortColumnIndex = -1;
	boolean ascending = true;

	public MetadataTableModel(GalleryData data) throws Exception
	{
		super(data);
		cols = visibleLabels.size();
		rows = n;
		renderer.hackBorders = false;
	}

	@Override
	public String getColumnName(int column)
	{
		return visibleLabels.get(column).getLabelName();
	}

	@Override
	public Class getColumnClass(int column)
	{
		ColumnInfo ci = visibleLabels.get(column);
		if (ci.render)
			return ImageItem.class;
		else if (ci.getLabel() == MDLabel.MDL_ENABLED)
			return Boolean.class;
		try
		{
			return MetaData.getLabelClass(ci.getLabel());
		}
		catch (Exception e)
		{
			e.printStackTrace();
			return null;
		}
	}

	@Override
	public int getIndex(int row, int col)
	{
		return row;
	}

	@Override
	public int[] getCoords(int index)
	{
		int[] coords = new int[2];
		coords[0] = index;
		coords[1] = 0;
		return coords;
	}

	@Override
	public Object getValueAt(int row, int column)
	{
		// DEBUG.printMessage(String.format("MetadataTable.getValueAt(%d, %d)",
		// row, column));
		try
		{
			ColumnInfo ci = visibleLabels.get(column);
			if (ci.render)
			{
				String key = getItemKey(row, ci.getLabel());
				ImageItem item;
				// If the element is on cache, just return it
				if (cache.containsKey(key))
					item = cache.get(key);
				else
				{
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
			switch (type)
			{
			case MetaData.LABEL_INT:
				int value = md.getValueInt(label, id);
				// treat special case of MDL_ENABLED
				if (label == MDLabel.MDL_ENABLED)
					return (value > 0);
				return value;
			case MetaData.LABEL_BOOL:
				return md.getValueBoolean(label, id);
			case MetaData.LABEL_DOUBLE:
				return md.getValueDouble(label, id);
			case MetaData.LABEL_SIZET:
				return md.getValueLong(label, id);
			case MetaData.LABEL_STRING:
				return md.getValueString(ci.getLabel(), data.ids[row]);
			case MetaData.LABEL_VECTOR_DOUBLE:
			case MetaData.LABEL_VECTOR_SIZET:
				return md.getValueString(ci.getLabel(), data.ids[row]);

			}
			return null;

		}
		catch (Exception e)
		{
			e.printStackTrace();
		}

		return null;
	}// function getValueAt

	@Override
	public void setValueAt(Object value, int row, int column)
	{

		ColumnInfo ci = visibleLabels.get(column);
		if (value == null || value.equals(""))
		{
			XmippDialog.showError(null, XmippMessage.getEmptyFieldMsg(ci.getLabelName()));
			return;//empty values are not allowed
		}
		try
		{
			long id = data.ids[row];
			setMdValueAt(value, row, column, ci, id);
			data.setMdChanges(true);
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}

	}// function setValueAt

	/**
	 * Helper function to set the value of a cell
	 * 
	 * @throws Exception
	 */
	protected void setMdValueAt(Object value, int row, int column, ColumnInfo ci, long id) throws Exception
	{
		if (!ci.render)
		{
			int label = ci.getLabel();
			int type = MetaData.getLabelType(label);
			MetaData md = data.md;
			switch (type)
			{
			case MetaData.LABEL_BOOL:
				md.setValueBoolean(label, (Boolean) value, id);
				break;
			case MetaData.LABEL_INT:
				if (label == MDLabel.MDL_ENABLED)
				{
					md.setEnabled((Boolean) value, id);
				}
				else
					md.setValueInt(label, ((Integer) value).intValue(), id);
				break;
			case MetaData.LABEL_DOUBLE:
				md.setValueDouble(label, ((Double) value).doubleValue(), id);
				break;
			case MetaData.LABEL_SIZET:
				md.setValueInt(label, ((Long) value).intValue(), id);
				break;
			case MetaData.LABEL_STRING:
			case MetaData.LABEL_VECTOR_DOUBLE:
			case MetaData.LABEL_VECTOR_SIZET:
				// TODO: Implement a better editor for vectors
				md.setValueString(label, value.toString(), id);
				break;
			}
			fireTableRowsUpdated(row, row);
		}
	}

	@Override
	public boolean isCellEditable(int row, int column)
	{
		ColumnInfo ci = visibleLabels.get(column);
		return ci.allowEdit && !ci.render;
	}

	// @Override
	// public String getImageFilenameAt(int row, int col){
	// ColumnInfo ci = visibleLabels.get(col);
	// return (ci.render && data.isImageFile(ci))
	// ? data.getValueFromCol(row, ci) : null;
	// }
	@Override
	public boolean handleDoubleClick(int row, int col)
	{
		try
		{
			ColumnInfo ci = visibleLabels.get(col);
			if (ci.allowRender && data.isImageFile(ci))
			{
				//new XmippImageWindow(data.window, new MdRowImageLoader(row, ci.getLabel()));
				ImagePlusLoader loader = new MdRowImageLoader(row, ci.getLabel());
				if(getNormalized())
					loader.setNormalize(normalize_min, normalize_max);
				ImagesWindowFactory.openXmippImageWindow(data.window, loader, loader.allowsPoll());
				return true;
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		return false;
	}// function handleDoubleClick

	@Override
	public void setColumns(int cols)
	{

	}

	@Override
	public boolean adjustColumn(int width)
	{
		return false;
	}

	@Override
	/** Whether to display the labels */
	public void setRenderImages(boolean value)
	{
		boolean changed = false;
		for (ColumnInfo ci : data.labels)
			if (ci.allowRender && ci.render != value)
			{
				ci.render = value;
				changed = true;
			}
		if (changed)
		{
			data.globalRender = value;
			calculateCellSize();
			fireTableDataChanged();
		}
	}

	@Override
	public String getLabel(int row, int col)
	{
		try
		{
			long objId = data.ids[row];
			return data.md.getValueString(visibleLabels.get(col).getLabel(), objId);
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		return null;
	}

	@Override
	protected void calculateCellSize()
	{
		// DEBUG.printMessage(String.format("MetadataTable:calculateSize"));
		if (data.globalRender)
		{
			super.calculateCellSize();
			// DEBUG.printMessage(String.format(
			// "MetadataTable:calculateSize w:%d, h:%d", cellDim.width,
			// cellDim.height));

		}
		else
		{
			int font_height;
			font_height = renderer.getFontMetrics(renderer.getFont()).getHeight();
			font_height += renderer.getIconTextGap(); // Adds the extra gap.
			cellDim.setSize(100, font_height + 5);
		}
	}

	@Override
	public void setupTable(JTable table)
	{
		table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		table.setDefaultRenderer(ImageItem.class, renderer);
		table.setDefaultRenderer(Double.class, new FloatRenderer());
		//MetadataTableHeader header = new MetadataTableHeader(columnModel);
		//JTableHeader header = new JTableHeader(table.getColumnModel());
		//table.setTableHeader(header);
		JTableHeader header = table.getTableHeader();
		header.setUpdateTableInRealTime(true);
		header.addMouseListener(new MetadataColumnListener(table));
		//header.setReorderingAllowed(true);
		updateTableSelection(table);
	}

	@Override
	public JTableHeader getTableHeaderModel()
	{
		return new MetadataTableHeader(columnModel);
	}

	/** Update the table selection according with data selection */
	@Override
	public void updateTableSelection(JTable table)
	{
		table.clearSelection();
		for (int i = 0; i < n; ++i)
			if (data.selection[i])
			{
				table.addRowSelectionInterval(i, i);
			}
	}

	@Override
	public boolean handleRightClick(int row, int col, XmippPopupMenuCreator xpopup)
	{
		xpopup.initItems();

		if (data.isFile(visibleLabels.get(col)))
		{
			xpopup.setItemVisible(XmippPopupMenuCreator.OPEN, true);
			if (!data.isImageFile(visibleLabels.get(col)))
				xpopup.setItemVisible(XmippPopupMenuCreator.OPEN_ASTEXT, true);
		}
		return true;
	}

	@Override
	/** Return the column model to be used with this table model */
	public GalleryColumnModel createColumnModel()
	{
		return new MetadataColumnModel();
	}

	public class MetadataColumnModel extends GalleryColumnModel
	{
		public MetadataColumnModel()
		{
			super(0);
		}

		@Override
		public void adjustColumnsWidth(JTable table)
		{
			try
			{
				if (visibleLabels.size() != getColumnCount())
					return;
				calculateCellSize();
				// String[] row = md.getRowValues(data.ids[0]);
				// table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
				int width = 0;
				TableCellRenderer rend;
				Component comp;
				boolean non_empty = data.md.size() > 0;

				for (int i = 0; i < visibleLabels.size(); ++i)
				{
					ColumnInfo col = visibleLabels.get(i);
					width = 0;
					// Calculate width of the cell
					if (col.render)
					{
						width = cellDim.width;
					}
					else if (non_empty)
					{
						// else {
						rend = table.getCellRenderer(0, i);
						comp = rend.getTableCellRendererComponent(table, getValueAt(0, i), false, false, 0, 0);
						width = comp.getPreferredSize().width + 10;
					}
					// Calculate width of the header
					TableColumn tc = getColumn(i);
					rend = tc.getHeaderRenderer();
					if (rend == null)
						rend = table.getTableHeader().getDefaultRenderer();
					Object value = tc.getHeaderValue();
					comp = rend.getTableCellRendererComponent(table, value, false, false, 0, i);
					// Take max width
					width = Math.max(width, comp.getPreferredSize().width);
					getColumn(i).setPreferredWidth(width + 10);
					// DEBUG.printMessage(String.format("col: %d, width: %d", i,
					// width));
				}
			}
			catch (Exception e)
			{
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}// function adjustColumnsWidth

	}// class MetadataColumnModel

	public class MetadataTableHeader extends JTableHeader
	{
		public MetadataTableHeader(TableColumnModel columnModel)
		{
			super(columnModel);
		}

		/** Show tooltips on columns header */
		public String getToolTipText(MouseEvent e)
		{
			java.awt.Point p = e.getPoint();
			int index = columnModel.getColumnIndexAtX(p.x);
			if (index > -1)
				return visibleLabels.get(index).comment;
			return null;
		}
	}

	public class MetadataColumnListener extends MouseAdapter
	{
		protected JTable table;

		public MetadataColumnListener(JTable t)
		{
			table = t;
		}


		public void mouseClicked(MouseEvent e)
		{
			TableColumnModel colModel = table.getColumnModel();
			//Get the clicked column index
			int columnModelIndex = colModel.getColumnIndexAtX(e.getX());
			//Take into account a possible reordering of columns
			int modelIndex = colModel.getColumn(columnModelIndex).getModelIndex();
			//Take into account possible invisible columns indexing
			modelIndex = data.getVisibleColumnIndex(modelIndex);
			if (modelIndex < 0)
				return;
			if (sortColumnIndex == modelIndex)
				ascending = !ascending;
			else
				sortColumnIndex = modelIndex;
			data.sortMd(sortColumnIndex, ascending);
			clearSelection();
			updateTableSelection(table);
			cache.clear();
			
			//fireTableDataChanged();
			
					
		}
	}

}// class MetadataTable
