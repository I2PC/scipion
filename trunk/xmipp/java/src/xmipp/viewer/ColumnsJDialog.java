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

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;

import xmipp.jni.MetaData;
import xmipp.utils.WindowUtil;

public class ColumnsJDialog extends JDialog implements ActionListener {
	private static final long serialVersionUID = 1L;
	
	private JFrameGallery parent;
	private JTable tableColumns;
	private JButton btnCancel;
	private JButton btnOk;
	private JButton btnUp;
	private JButton btnDown;
	private ColumnsTableModel model;
	// This will be used for check for results from the dialog
	private ArrayList<ColumnInfo> rows;
	boolean fireEvent = true;

	public ColumnsJDialog(JFrameGallery parent, boolean modal) {
		super(parent, modal);
		this.parent = parent;
		initComponents();
	}// constructor ColumnsJDialog

	public ArrayList<ColumnInfo> getColumnsResult() {
		return rows;
	}

	private void initComponents() {
		setResizable(false);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Columns");
		setLayout(new GridBagLayout());
		GridBagConstraints gbc = new GridBagConstraints();
		gbc.insets = new Insets(10, 10, 10, 10);
		gbc.anchor = GridBagConstraints.WEST;

		JPanel groupstbpn = new JPanel();
		JScrollPane sp = new JScrollPane();
		groupstbpn.setBorder(BorderFactory
				.createTitledBorder("Column properties"));
		groupstbpn.add(sp);
		sp.setOpaque(true);
		model = new ColumnsTableModel(parent.getData().labels);
		tableColumns = new JTable(model);
		tableColumns
				.setPreferredScrollableViewportSize(new Dimension(350, 200));
		tableColumns
				.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		sp.setViewportView(tableColumns);
		add(groupstbpn, WindowUtil.getConstraints(gbc, 0, 0, 3, 3));
		btnCancel = WindowUtil.getTextButton("Cancel", this);
		btnOk = WindowUtil.getTextButton("Ok", this);
		JPanel panel = new JPanel();
		panel.setLayout(new GridBagLayout());
		btnUp = WindowUtil.getIconButton("up.gif", this);

		panel.add(btnUp, WindowUtil.getConstraints(gbc, 0, 0));
		btnDown = WindowUtil.getIconButton("down.gif", this);
		panel.add(btnDown, WindowUtil.getConstraints(gbc, 0, 1));
		add(panel, WindowUtil.getConstraints(gbc, 3, 2));
		add(btnCancel, WindowUtil.getConstraints(gbc, 3, 3));
		gbc.anchor = GridBagConstraints.LINE_END;
		add(btnOk, WindowUtil.getConstraints(gbc, 2, 3));
		getRootPane().setDefaultButton(btnOk);

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

		pack();
		WindowUtil.centerWindows(this, parent);
		setVisible(true);
	}// function initComponents

	private void enableUpDown(boolean value) {
		btnUp.setEnabled(value);
		btnDown.setEnabled(value);
	}// function enableUpDown

	private void close() {
		setVisible(false);
		dispose();
	}// function close

	// move the selection on the table, -1 up, 0 down
	private void moveSelection(int deltha) {
		int pos = tableColumns.getSelectedRow();
		ColumnInfo ci = rows.remove(pos);
		pos += deltha;
		rows.add(pos, ci);
		model.fireTableDataChanged();
		tableColumns.setRowSelectionInterval(pos, pos);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		JButton btn = (JButton) e.getSource();

		if (btn == btnUp && tableColumns.getSelectedRow() > 0)
			moveSelection(-1);
		else if (btn == btnDown && tableColumns.getSelectedRow() < rows.size() - 1)
			moveSelection(1);
		else if (btn == btnOk)
			close();
		else if (btn == btnCancel) {
			rows = null;
			close();
		}
	}// function actionPerformed

	class ColumnsTableModel extends AbstractTableModel {
		private static final long serialVersionUID = 1L;

		private String[] columns = { "Label", "Visible", "Render" };

		public ColumnsTableModel(int[] labels) {
			rows = new ArrayList<ColumnInfo>(labels.length);
			for (int i = 0; i < labels.length; ++i)
				rows.add(new ColumnInfo(labels[i]));
		}

		public ColumnsTableModel(ArrayList<ColumnInfo> labelsInfo) {
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
				if (column == 0
						|| (column == 2 && !MetaData.isImage(rows.get(row)
								.getLabel())))
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
			case 1:
				col.visible = (Boolean) value;
				break;
			case 2:
				col.render = (Boolean) value;
				break;
			}

		}

		@Override
		public Object getValueAt(int row, int column) {
			ColumnInfo col = rows.get(row);
			switch (column) {
			case 0:
				return col.getLabelName();
			case 1:
				return col.visible;
			case 2:
				return col.render;
			}
			return null;
		}

	}// class ColumnsTableModel

}// class ColumnsJDialog
