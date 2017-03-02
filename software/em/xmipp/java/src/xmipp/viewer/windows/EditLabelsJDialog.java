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

import xmipp.utils.XmippDialog;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.models.ImageGalleryTableModel;

public class EditLabelsJDialog extends XmippDialog {
	private static final long serialVersionUID = 1L;
	private JTable tableColumns;
	private JButton btnAdd, btnDelete, btnFill, btnOpen;
	private List<ColumnInfo> rows;
	private ColumnsTableModel model;
	// This will be used for check for results from the dialog
	boolean fireEvent = true;
	GridBagConstraints gbc = new GridBagConstraints();
	ImageGalleryTableModel gallery;
	JPanel panelButtons;

	public EditLabelsJDialog(GalleryJFrame parent) {
		super(parent, "Edit labels", true);
		this.rows = parent.getData().getLabelsInfo();
		this.gallery = parent.gallery;
		btnOkText = "Close";
		btnCancelDisplay = false;
		initComponents();
		enableDelete(false);
		enableFill(false);
	}// constructor ColumnsJDialog

	private JButton createButton(String icon, String tip) {
		JButton btn = XmippWindowUtil.getIconButton(icon, this);
		btn.setToolTipText(tip);
		btn.setFocusable(false);
		panelButtons.add(btn);
		return btn;
	}

	protected void createToolbarButtons() {
		panelButtons = new JPanel();
		btnAdd = createButton("add.gif", "Add new label");
		btnDelete = createButton("delete.gif", "Delete label");
		btnFill = createButton("fill.png", "Fill label");
	}

	@Override
	protected void createContent(JPanel panel) {
		setResizable(false);
		panel.setLayout(new GridBagLayout());
		gbc.anchor = GridBagConstraints.EAST;

		JPanel groupstbpn = new JPanel();
		JScrollPane sp = new JScrollPane();
		groupstbpn.setBorder(BorderFactory.createTitledBorder("Labels"));
		groupstbpn.add(sp);
		sp.setOpaque(true);
		model = new ColumnsTableModel();
		tableColumns = new JTable(model);
		tableColumns
				.setPreferredScrollableViewportSize(new Dimension(350, 200));
		sp.setViewportView(tableColumns);
		panel.add(groupstbpn, XmippWindowUtil.getConstraints(gbc, 0, 1, 2));
		createToolbarButtons();
		panel.add(panelButtons, XmippWindowUtil.getConstraints(gbc, 1, 0));

		// listen to selection changes (only one row selected)
		// tableColumns.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		tableColumns.getSelectionModel().addListSelectionListener(
				new ListSelectionListener() {
					@Override
					public void valueChanged(ListSelectionEvent e) {
					    int n = getSelectedRows().length;
						enableDelete(n > 0);
						enableFill(n == 1);
					}
				});
	}// function initComponents

	private void enableDelete(boolean value) {
		btnDelete.setEnabled(value);
	}

	private void enableFill(boolean value) {
	    btnFill.setEnabled(value);
	}

	private void rowsChanged() {
		this.rows = ((GalleryJFrame) parent).getData().getLabelsInfo();
		model.fireTableDataChanged();
	}

	private void showAddFillDialog(AddFillLabelsJDialog dlg) throws Exception {

		if (dlg.showDialog()) {
			((GalleryJFrame) parent).fillLabel(dlg.getLabel(),
					dlg.getFillMode(), dlg.getValues());
			rowsChanged();
		}

	}

	@Override
	public void handleActionPerformed(ActionEvent evt) {
		try {
			JButton btn = (JButton) evt.getSource();
			if (btn == btnAdd) {
				showAddFillDialog(new AddFillLabelsJDialog(
						(GalleryJFrame) parent, rows));
			} else if (btn == btnDelete) {
			    System.out.println("Here testing java...");
				if (showWarning("Are you sure to remove this label?")) {
				    int [] rows = getSelectedRows();
				    System.out.println("Rows: " + Integer.toString(rows.length));
				    for (int row: rows) {
				        ((GalleryJFrame) parent).removeLabel(rowToLabel(row));
				    }
					rowsChanged();
				}
			} else if (btn == btnFill) {
			    int [] rows = getSelectedRows();
				showAddFillDialog(new AddFillLabelsJDialog(
						(GalleryJFrame) parent, rowToLabel(rows[0])));
			}
		} catch (Exception e) {
			showException(e);
		}
	}// function actionPerformed

	public int[] getSelectedRows() {
		return tableColumns.getSelectedRows();
	}

    // Get the label corresponding to a given row
    public int rowToLabel(int row) {
        return rows.get(row).label;
    }

	class ColumnsTableModel extends AbstractTableModel {
		private static final long serialVersionUID = 1L;

		private String[] columns = { "Label", "Type" };

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
			return false;
		}

		@Override
		public void setValueAt(Object value, int row, int column) {
		}

		@Override
		public Object getValueAt(int row, int column) {
			ColumnInfo col = rows.get(row);
			switch (column) {
			case 0:
				return col.labelName;
			case 1:
				return col.getLabelTypeString();
			}
			return null;
		}

	}// class ColumnsTableModel

}// class ColumnsJDialog
