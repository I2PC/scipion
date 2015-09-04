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


import javax.swing.JTable;

import xmipp.jni.MetaData;
import xmipp.utils.Params;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.viewer.windows.GalleryJFrame;

/** This is a data model designed for Row metadatas */
@SuppressWarnings("serial")
public class MetadataRowTableModel extends MetadataTableModel {
	protected long id;
	
	/** Constructor using a metadata row */
	public MetadataRowTableModel(GalleryJFrame window, MetaData md) throws Exception {
		this(new GalleryData(window, new Params(), md));
	}
	
	public MetadataRowTableModel(GalleryData data) throws Exception {
		super(data, null);
		cols = 1;
		rows = visibleLabels.size();
		id = data.md.firstObject();
		// TODO Auto-generated constructor stub
                selection = new boolean[rows];
	}
	
	@Override
	public String getColumnName(int column) {
		return "Value";
	}

	@Override
	public Class getColumnClass(int column) {
		return String.class;
	}
	
	@Override
	public Object getValueAt(int row, int column){
		return super.getValueAt(0, row);
	}
	
	@Override
	public void setValueAt(Object value, int row, int column) {
		data.setValueToCol(0, data.getColumnInfo(row), value.toString());
	}// function setValueAt
	
	@Override
	/** Return the column model to be used with this table model */
	public GalleryColumnModel createColumnModel() {
		return new GalleryColumnModel(cellDim.width);
	}//function createColumnModel

	@Override
	public void setupTable(JTable table) {
		super.setupTable(table);
		table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
	}
	
	@Override
	public boolean handleRightClick(int row, int col, XmippPopupMenuCreator xpopup) {
		boolean result =  super.handleRightClick(0, row, xpopup);
                xpopup.setItemVisible(XmippPopupMenuCreator.SELECT, false);
                xpopup.setItemVisible(XmippPopupMenuCreator.ENABLED, false);
                xpopup.setItemVisible(XmippPopupMenuCreator.DISABLED, false);
                return result;
	}
        
       @Override
	public int getIndex(int row, int col) {

                    return 0;

	}
       
    public ColumnInfo getColumn(int row, int col)
   	{
   		return visibleLabels.get(row);
   	}
        
        /** Set the selection state of an element give row and col */
        @Override
	public void touchItem(int row, int col) {
                    setSelected(row, !isSelected(row));
                    adjustWidth = false;
                    //fireTableCellUpdated(row, col);//item was clicked
	}
        
    public boolean isSelected(int row, int col) {
       
        return isSelected(row);
    }
    
    public void setSelected(int row, int col, boolean b) {
        setSelected(row, b);
    }

	//** Select a range of elements given the coordinates */
        @Override
	public void selectRange(int first_row, int first_col, int last_row, int last_col) {
		int min = Math.min(first_row, last_row);
                int max = Math.max(first_row, last_row);
		for (int i = min; i <= max; i++)
            if(!isSelected(i))
            {
				setSelected(i, true);
	        	fireTableCellUpdated(i, 0);
            }
	}
        
          @Override
    public boolean showLabels() {
        return false;
    }
}//class MetadataRow
