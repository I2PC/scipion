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

package xmipp.viewer.windows;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;

import xmipp.ij.commons.XmippApplication;
import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.models.ColumnInfo;

public class ColumnsJDialog extends XmippDialog {
	private static final long serialVersionUID = 1L;
	private JTable tableColumns;
	private JButton btnUp;
	private JButton btnDown;
	private ColumnsTableModel model;
	// This will be used for check for results from the dialog
	private List<ColumnInfo> rows;
	boolean fireEvent = true;

	public ColumnsJDialog(GalleryJFrame parent) {
		super(parent, "Columns", true);
		initComponents();
	}// constructor ColumnsJDialog

	public List<ColumnInfo> getColumnsResult() {
		return rows;
	}

	@Override
	protected void createContent(JPanel panel){
		setResizable(false);
		panel.setLayout(new GridBagLayout());
		GridBagConstraints gbc = new GridBagConstraints();
		gbc.insets = new Insets(10, 10, 10, 10);
		gbc.anchor = GridBagConstraints.WEST;

		JPanel groupstbpn = new JPanel();
		JScrollPane sp = new JScrollPane();
		groupstbpn.setBorder(BorderFactory
				.createTitledBorder("Column properties"));
		groupstbpn.add(sp);
		sp.setOpaque(true);
		model = new ColumnsTableModel(((GalleryJFrame)parent).getData().getLabelsInfo());
		tableColumns = new JTable(model);
		tableColumns
				.setPreferredScrollableViewportSize(new Dimension(350, 200));
		sp.setViewportView(tableColumns);
		panel.add(groupstbpn, XmippWindowUtil.getConstraints(gbc, 0, 0));

		JPanel panelUpDown = new JPanel();
		panelUpDown.setLayout(new GridBagLayout());
		gbc.insets = new Insets(0, 0, 5, 5);
		btnUp = XmippWindowUtil.getIconButton("up.gif", this);
		panelUpDown.add(btnUp, XmippWindowUtil.getConstraints(gbc, 0, 0));
		btnDown = XmippWindowUtil.getIconButton("down.gif", this);
		panelUpDown.add(btnDown, XmippWindowUtil.getConstraints(gbc, 0, 1));
		panel.add(panelUpDown, XmippWindowUtil.getConstraints(gbc, 1, 0));
		// this buttons will be enabled after selection
		enableUpDown(false);
		// listen to selection changes (only one row selected)
		tableColumns.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		tableColumns.getSelectionModel().addListSelectionListener(
				new ListSelectionListener() {
					@Override
					public void valueChanged(ListSelectionEvent e) {
						enableUpDown(true);
					}
				});
		formatTable();
	}// function initComponents
	
	protected void formatTable() {
        tableColumns.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        tableColumns.getColumnModel().getColumn(0).setPreferredWidth(200);
        tableColumns.getColumnModel().getColumn(1).setPreferredWidth(50);
        tableColumns.getColumnModel().getColumn(2).setPreferredWidth(50);
        tableColumns.getColumnModel().getColumn(3).setPreferredWidth(50);
        
    }

	protected void enableUpDown(boolean value) {
		btnUp.setEnabled(value);
		btnDown.setEnabled(value);
	}// function enableUpDown

	// move the selection on the table, -1 up, 0 down
	protected void moveSelection(int deltha) {
		int pos = tableColumns.getSelectedRow();
		ColumnInfo ci = rows.remove(pos);
		pos += deltha;
		rows.add(pos, ci);
		model.fireTableDataChanged();
		tableColumns.setRowSelectionInterval(pos, pos);
	}

	@Override
	public void handleActionPerformed(ActionEvent evt){		
		JButton btn = (JButton) evt.getSource();
		
		if (btn == btnUp && tableColumns.getSelectedRow() > 0)
			moveSelection(-1);
		else if (btn == btnDown && tableColumns.getSelectedRow() < rows.size() - 1)
			moveSelection(1);
	}// function actionPerformed

	class ColumnsTableModel extends AbstractTableModel {
		private static final long serialVersionUID = 1L;

		private String[] columns = { "Label", "Visible", "Render", "Edit" };

		public ColumnsTableModel(int[] labels) {
			rows = new ArrayList<ColumnInfo>(labels.length);
			for (int i = 0; i < labels.length; ++i)
				rows.add(new ColumnInfo(labels[i]));
		}

		public ColumnsTableModel(List<ColumnInfo> labelsInfo) {
			int n = labelsInfo.size();
			rows = new ArrayList<ColumnInfo>(n);
			for (int i = 0; i < n; ++i)
				rows.add(labelsInfo.get(i).clone());
		}

		@Override
		public Class getColumnClass(int column) {
			return getValueAt(0, column).getClass();
		}

		@Override
		public String getColumnName(int columnIndex) {
			return columns[columnIndex];
		}

		@Override
		public int getColumnCount() {
			return columns.length;
		}

		@Override
		public int getRowCount() {
			return rows.size();
			// return frame.getParticlePicker().getFamilies().size();
		}

		@Override
		public boolean isCellEditable(int row, int column) {
			try {
				if (column == 2 && !rows.get(row).allowRender)
					return false;
				GalleryJFrame frame = (GalleryJFrame)parent;
				
				if (column == 3 && frame.data.isScipionInstance())
					return false;
			} catch (Exception e) {
				e.printStackTrace();
			}
			return true;
		}

		@Override
		public void setValueAt(Object value, int row, int column) {
			ColumnInfo col = rows.get(row);
			switch (column) {
			case 0:
				col.labelName = ((String)value);
				break;
			case 1:
				col.visible = (Boolean) value;
				break;
			case 2:
				col.render = (Boolean) value;
				break;
			case 3:
				col.allowEdit = (Boolean) value;
				break;
			}

		}

		@Override
		public Object getValueAt(int row, int column) {
			ColumnInfo col = rows.get(row);
			switch (column) {
			case 0:
				return col.labelName;
			case 1:
				return col.visible;
			case 2:
				return col.render;
			case 3:
				return col.allowEdit;
			}
			return null;
		}

	}// class ColumnsTableModel

}// class ColumnsJDialog
