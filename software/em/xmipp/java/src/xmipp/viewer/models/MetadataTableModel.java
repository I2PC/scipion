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
import java.awt.Point;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.*;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;

import xmipp.ij.commons.XmippApplication;
import xmipp.ij.commons.XmippUtil;
import xmipp.jni.Filename;
import xmipp.jni.MetaData;
import xmipp.utils.Params;
import xmipp.utils.ScipionParams;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.FloatRenderer;
import xmipp.viewer.windows.ImagesWindowFactory;

public class MetadataTableModel extends MetadataGalleryTableModel {
	private static final long serialVersionUID = 1L;
    final static int NO_COLUMN_INDEX = -1;
	int sortColumnIndex = NO_COLUMN_INDEX;
	boolean ascending = true;
        

	public MetadataTableModel(GalleryData data, boolean[] selection) throws Exception {
		super(data, selection);
		cols = visibleLabels.size();
		rows = n;
		renderer.hackBorders = false;
                
	}

	@Override
	public String getColumnName(int column) {
		return visibleLabels.get(column).labelName;
	}

	@Override
	public Class getColumnClass(int column) {
		ColumnInfo ci = visibleLabels.get(column);
		Class c = null;
		try {
			if (ci.render)
				c = ImageItem.class;
			else if (ci.isEnable())
				c = Boolean.class;//This way a JCheckBox is rendered
			else
				c = MetaData.getClassForType(ci.type);
                        
		} catch (Exception e) {
			e.printStackTrace();
		}
		return c;
	}

	
	/**
	 * Returns metadata value with java type
	 */
	@Override
	public Object getValueAt(int row, int column) {
		// DEBUG.printMessage(String.format("MetadataTable.getValueAt(%d, %d)",
		// row, column));
		try {
                        
			ColumnInfo ci = visibleLabels.get(column);
			if (ci.render) {
				
				String key = getItemKey(row, ci.label);
				ImageItem item;
				// If the element is on cache, just return it
				if (cache.containsKey(key))
					item = cache.get(key);
				else {
					// If not, create the item and store it for future
					item = createImageItem(row, ci);
					cache.put(key, item);
				}
				setupItem(item);
				return item;
			}
			int label = ci.label;
			long id = data.ids[row];
                        
			int type = ci.type;
			MetaData md = data.md;
			switch (type) {
				case MetaData.LABEL_INT:
	                            
					int value = md.getValueInt(label, id);
					// treat special case of MDL_ENABLED
					if (ci.isEnable())
						return (value > 0);
					return value;
				case MetaData.LABEL_BOOL:
					return md.getValueBoolean(label, id);
				case MetaData.LABEL_DOUBLE:
					return md.getValueDouble(label, id);
				case MetaData.LABEL_SIZET:
					return md.getValueLong(label, id);
				case MetaData.LABEL_STRING:
                    String str = md.getValueString(label, data.ids[row]);
                    if (ci.labelName.contains("_transform._matrix"))
                        return String.format("<html>%s</html>", XmippUtil.formatNumbers(str).replace("],", "]<br>"));
	
					return str;
				case MetaData.LABEL_VECTOR_DOUBLE:
				case MetaData.LABEL_VECTOR_SIZET:
					return md.getValueString(label, data.ids[row]);

			}
			return null;

		} catch (Exception e) {
                    e.printStackTrace();
                    
		}
                return null;
	}// function getValueAt
        
       

	@Override
	public void setValueAt(Object value, int row, int column) {

		ColumnInfo ci = visibleLabels.get(column);
		if (value == null || value.equals("")) {
			XmippDialog.showError(null,
					XmippMessage.getEmptyFieldMsg(ci.labelName));
			return;// empty values are not allowed
		}
		try {
			long id = data.ids[row];
			setMdValueAt(value, row, column, ci, id);
			data.setMdChanges(true);
		} catch (Exception e) {
			e.printStackTrace();
		}

	}// function setValueAt

	/**
	 * Helper function to set the value of a cell
	 * 
	 * @throws Exception
	 */
	protected void setMdValueAt(Object value, int row, int column,
			ColumnInfo ci, long id) throws Exception {
		if (!ci.render) {
			int label = ci.label;
			int type = ci.type;
			MetaData md = data.md;
			switch (type) {
			case MetaData.LABEL_BOOL:
				md.setValueBoolean(label, (Boolean) value, id);
				break;
			case MetaData.LABEL_INT:
				if (ci.isEnable()) {
					md.setEnabled((Boolean) value, id);
				} else
					md.setValueInt(label, ((Integer) value).intValue(), id);
				break;
			case MetaData.LABEL_DOUBLE:
				md.setValueDouble(label, ((Double) value).doubleValue(), id);
				break;
			case MetaData.LABEL_SIZET:
				md.setValueLong(label, ((Long) value).longValue(), id);
				break;
			case MetaData.LABEL_STRING:
			case MetaData.LABEL_VECTOR_DOUBLE:
			case MetaData.LABEL_VECTOR_SIZET:
				// TODO: Implement a better editor for vectors
				md.setValueString(label, value.toString(), id);
				break;
			}
			adjustWidth = false;
			fireTableRowsUpdated(row, row);
		}
	}

	@Override
	public boolean isCellEditable(int row, int column) {
        ColumnInfo ci = visibleLabels.get(column);
        if(ci.isEnable())
            return ci.allowEdit;//maybe on metadatas sometimes is disabled
        if(data.isScipionInstance() || data.isChimeraClient())
                return false;
		
		return ci.allowEdit && !ci.render;
	}

	// @Override
	// public String getImageFilenameAt(int row, int col){
	// ColumnInfo ci = visibleLabels.get(col);
	// return (ci.render && data.isImageFile(ci))
	// ? data.getValueFromCol(row, ci) : null;
	// }
	@Override
	public boolean handleDoubleClick(int row, int col) {
		try {
			
			ColumnInfo ci = visibleLabels.get(col);
			if (ci.allowRender && data.isImageFile(ci)) {
                int index = getIndex(row, col);
                String file = getImageFilename(index, ci.label);
                if(Filename.isVolume(file))
                {
                	Params params = XmippApplication.isScipion()? ((ScipionParams)data.parameters).getScipionParams(): new Params();
                    ImagesWindowFactory.openMetadata(file, params, Params.OPENING_MODE_GALLERY);
                }
                else
                    openXmippImageWindow(index, ci);
				return true;
			}
            if(data.isChimeraClient())//
            {
                String cmd = data.getChimeraProjectionCmd(row);
                XmippWindowUtil.runCommand(cmd, data.parameters.getChimeraPort());
            }
		} catch (Exception e) {
            if(e.getMessage() == null)
                e.printStackTrace();
            else
                XmippDialog.showError(null, e.getMessage());
			
		}
		return false;
	}// function handleDoubleClick

	@Override
	public void setColumns(int cols) {

	}

	@Override
	public boolean adjustColumn(int width) {
		return false;
	}

	/** Whether to display the labels */
	public void setRenderImages(boolean value) {
		boolean changed = false;
		for (ColumnInfo ci : data.labels)
			if (ci.allowRender && ci.render != value) {
				ci.render = value;
				changed = true;
				if(data.ciFirstRender == null && value)
					data.ciFirstRender = ci;
			}
		try
		{
			if (changed) {
				data.renderImages = value;
				loadDimension();
				calculateCellSize();
				fireTableDataChanged();
			}
			
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	@Override
        public boolean showLabels() {
            return data.renderImages;
        }

	@Override
	protected void calculateCellSize() {
		//System.out.println(String.format("MetadataTable:calculateSize renderLabel %s hasRenderLabel %s", data.renderImages, data.hasRenderLabel()));
		if (data.renderImages && data.hasRenderLabel()) {
			super.calculateCellSize();
			//System.out.println(String.format("MetadataTable:calculateSize w:%d, h:%d", cellDim.width, cellDim.height));

		} else {
			int font_height;
			font_height = renderer.getFontMetrics(renderer.getFont())
					.getHeight();
			font_height += renderer.getIconTextGap(); // Adds the extra gap.
			cellDim.setSize(200, font_height + 5);
		}
	}

	@Override
	public void setupTable(JTable table) {
		table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		table.setDefaultRenderer(ImageItem.class, renderer);
		table.setDefaultRenderer(Double.class, new FloatRenderer());
		// MetadataTableHeader header = new MetadataTableHeader(columnModel);
		// JTableHeader header = new JTableHeader(table.getColumnModel());
		// table.setTableHeader(header);
		JTableHeader header = table.getTableHeader();
		header.setUpdateTableInRealTime(true);
		header.addMouseListener(new MetadataColumnListener(table));
		// header.setReorderingAllowed(true);
		updateTableSelection(table);
	}

	@Override
	public JTableHeader getTableHeaderModel() {
		return new MetadataTableHeader(columnModel);
	}

	/** Update the table selection according with data selection */
	@Override
	public void updateTableSelection(JTable table) {
		table.clearSelection();
		for (int i = 0; i < selection.length; ++i)
			if (selection[i]) {
				table.addRowSelectionInterval(i, i);
			}
	}

	@Override
	public boolean handleRightClick(int row, int col,
			XmippPopupMenuCreator xpopup) {
             if (isBusy(row, col))
                    return false;
		xpopup.initItems();
                
		if (data.isFile(visibleLabels.get(col))) {
                        
			xpopup.setItemVisible(XmippPopupMenuCreator.OPEN, true);
			if (!data.isImageFile(visibleLabels.get(col)))
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
				if (visibleLabels.size() != getColumnCount())
					return;
				
				//DEBUG.printStackTrace();
				
				calculateCellSize();
				// String[] row = md.getRowValues(data.ids[0]);
				// table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
				int width = 0;
				TableCellRenderer rend;
				Component comp;
				boolean non_empty = data.md.size() > 0;
				TableColumn tc;
				ColumnInfo col;
				for (int i = 0; i < visibleLabels.size(); i++) {
					 col = visibleLabels.get(i);
					tc = getColumn(i);
					width = 0;
					// Calculate width of the cell
					if (col.render) {
						width = cellDim.width;
					} else if (non_empty) {
						// else {
						
						rend = table.getCellRenderer(0, i);
						if(rend != null)
						{
							Object value = table.getValueAt(0, i);
							//If columns are reordered Renderer may be related to column according to visual order but value is related to real order
							comp = rend.getTableCellRendererComponent(table, value, false, false, 0, i);
							width = comp.getPreferredSize().width + 10;
						}
					}
					// Calculate width of the header
					rend = tc.getHeaderRenderer();
					if (rend == null)
						rend = table.getTableHeader().getDefaultRenderer();
					Object value = tc.getHeaderValue();
					comp = rend.getTableCellRendererComponent(table, value,
							false, false, 0, i);
					// Take max width
					width = Math.max(width, comp.getPreferredSize().width);
					getColumn(i).setPreferredWidth(width + 10);
					// DEBUG.printMessage(String.format("col: %d, width: %d", i,
					// width));
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}// function adjustColumnsWidth

	}// class MetadataColumnModel

	public class MetadataTableHeader extends JTableHeader {
		public MetadataTableHeader(TableColumnModel columnModel) {
			super(columnModel);
		}

		/** Show tooltips on columns header */
		public String getToolTipText(MouseEvent e) {
			java.awt.Point p = e.getPoint();
			int index = columnModel.getColumnIndexAtX(p.x);
			if (index > -1){
                int col =  table.convertColumnIndexToModel(index);
                return visibleLabels.get(col).comment;
            }

			return null;
		}
	}

	public class MetadataColumnListener extends MouseAdapter {
        public static final String ASCENDING_PREFIX = "\u25B2 ";
        public static final String DESCENDING_PREFIX = "\u25BC ";
        protected JTable table;

		public MetadataColumnListener(JTable t) {
			table = t;
		}

		public void mouseClicked(MouseEvent e) {
                        
			TableColumnModel colModel = table.getColumnModel();
			// Get the clicked column index
			int columnModelIndex = colModel.getColumnIndexAtX(e.getX());

            // Get the column
            final TableColumn column = colModel.getColumn(columnModelIndex);

            // Take into account a possible reordering of columns
			int modelIndex = column.getModelIndex();

			// Take into account possible invisible columns indexing
			modelIndex = data.getVisibleColumnIndex(modelIndex);
			if (modelIndex < 0)
				return;
			if (sortColumnIndex == modelIndex)
				ascending = !ascending;
			else

                // Remove previous sorting icon
                removePreviousSortingIcon();
                sortColumnIndex = modelIndex;


            final String columnName = data.labels.get(sortColumnIndex).labelName;

            column.setHeaderValue( "\u231B " + columnName);
            Runnable sort = new Runnable() {
                @Override
                public void run() {

                    data.sortMd(data.labels.get(sortColumnIndex), ascending);
                    column.setHeaderValue((ascending? ASCENDING_PREFIX : DESCENDING_PREFIX)+columnName);
                    table.getTableHeader().repaint();
                    clearSelection();
                    updateTableSelection(table);
                    cache.clear();
                    table.repaint();

                }
            };

            SwingUtilities.invokeLater(sort);

		}

        protected void removePreviousSortingIcon() {

            // If there is a previous column sorted
            if (sortColumnIndex != NO_COLUMN_INDEX){
                // Get the column model
                TableColumnModel colModel = table.getColumnModel();

                // For each column
                for (int i=0; i<colModel.getColumnCount();i++){
                    TableColumn col = colModel.getColumn(i);
                    String header = col.getHeaderValue().toString();

                    // If it has the ascending prefix
                    if (header.startsWith(ASCENDING_PREFIX) ){
                        // remove it
                        col.setHeaderValue(header.replaceFirst(ASCENDING_PREFIX, ""));
                        return;
                    } else if (header.startsWith(DESCENDING_PREFIX)) {
                        // remove it
                        col.setHeaderValue(header.replaceFirst(DESCENDING_PREFIX, ""));
                        return;
                    }
                }
            }
        }
    }
        
	@Override
	public int getIndex(int row, int col) {
        return row;
	}
    
    public ColumnInfo getColumn(int row, int col)
	{
		return visibleLabels.get(col);
	}
    
    @Override
	public Point getCoords(int index) {
	
		Point p = new Point();
		p.x = 0;
		p.y = index;
		return p;
	}
	

        
        
}// class MetadataTable
