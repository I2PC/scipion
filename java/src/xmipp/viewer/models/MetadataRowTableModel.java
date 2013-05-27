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

import java.awt.Window;

import javax.swing.JTable;

import xmipp.ij.commons.XmippImageWindow;
import xmipp.jni.MetaData;
import xmipp.utils.Param;
import xmipp.utils.XmippPopupMenuCreator;

/** This is a data model designed for Row metadatas */
@SuppressWarnings("serial")
public class MetadataRowTableModel extends MetadataTableModel {
	protected long id;
	
	/** Constructor using a metadata row */
	public MetadataRowTableModel(Window window, MetaData md) throws Exception {
		this(new GalleryData(window, null, new Param(), md));
	}
	
	public MetadataRowTableModel(GalleryData data) throws Exception {
		super(data);
		cols = 1;
		rows = visibleLabels.size();
		id = data.md.firstObject();
		// TODO Auto-generated constructor stub
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
	public boolean handleRightClick(int row, int col,
			XmippPopupMenuCreator xpopup) {
		return super.handleRightClick(0, row, xpopup);
	}
	
}//class MetadataRow
