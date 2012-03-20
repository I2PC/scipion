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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.Icon;
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
import xmipp.utils.ColorEditor;
import xmipp.utils.ColorIcon;
import xmipp.utils.ColorRenderer;
import xmipp.utils.XmippWindowUtil;
import xmipp.utils.XmippDialog;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.models.ImageGallery;

public class ClassesJDialog extends XmippDialog {
	private static final long serialVersionUID = 1L;
	private JTable tableClasses;
	private JButton btnAdd, btnDelete, btnSave, btnOpen;
	private ClassesTableModel model;
	// This will be used for check for results from the dialog
	private ArrayList<ClassInfo> classes;
	boolean fireEvent = true;
	GridBagConstraints gbc = new GridBagConstraints();
	ImageGallery gallery;
	JPanel panelButtons;
	
	public ClassesJDialog(JFrameGallery parent, ArrayList<ClassInfo> classes) {
		super(parent, "Columns", true);
		this.classes = classes;
		this.gallery = parent.gallery;
		this.btnOkText = "Select";
		disposeOnClose = false;
		initComponents();
		enableDelete(false);
		
	}// constructor ColumnsJDialog

	private JButton createButton(String icon){
		JButton btn = XmippWindowUtil.getIconButton(icon, this);
		btn.setFocusable(false);
		panelButtons.add(btn);
		return btn;
	}
	
	protected void createToolbarButtons(){
		panelButtons = new JPanel();
		btnAdd = createButton("add.gif");
		btnDelete = createButton("delete.gif");
		btnSave = createButton("save.gif");
	}
	
	@Override
	protected void createContent(JPanel panel) {
		setResizable(false);
		panel.setLayout(new GridBagLayout());
		gbc.anchor = GridBagConstraints.EAST;

		JPanel groupstbpn = new JPanel();
		JScrollPane sp = new JScrollPane();
		groupstbpn.setBorder(BorderFactory.createTitledBorder("Classes"));
		groupstbpn.add(sp);
		sp.setOpaque(true);
		model = new ClassesTableModel();
		tableClasses = new JTable(model);
		tableClasses
				.setPreferredScrollableViewportSize(new Dimension(350, 200));
		sp.setViewportView(tableClasses);
		tableClasses.setDefaultRenderer(Color.class, new ColorRenderer());
		tableClasses.setDefaultEditor(Color.class, new ColorEditor());
		panel.add(groupstbpn, XmippWindowUtil.getConstraints(gbc, 0, 1, 2));
		createToolbarButtons();
		panel.add(panelButtons, XmippWindowUtil.getConstraints(gbc, 1, 0));
		
		// listen to selection changes (only one row selected)
		tableClasses.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		tableClasses.getSelectionModel().addListSelectionListener(
				new ListSelectionListener() {
					@Override
					public void valueChanged(ListSelectionEvent e) {
						enableDelete(getSelectedClass() >= 0);
					}
				});
	}// function initComponents
	
	/** Override this method to update the class and image count */
	@Override 
	public boolean showDialog(){
		gallery.data.updateClassesInfo();
		return super.showDialog();
	}

	private void enableDelete(boolean value) {
		//btnAdd.setEnabled(value);
		btnDelete.setEnabled(value);
		btnOk.setEnabled(value);
	}// function enableUpDown

	@Override
	public void handleActionPerformed(ActionEvent evt) {
		JButton btn = (JButton) evt.getSource();
		if (btn == btnAdd){
			model.addNewRow();
			int row = model.getRowCount() - 1;
			tableClasses.editCellAt(row, 0);
			tableClasses.setRowSelectionInterval(row, row);
		}
		else if (btn == btnDelete){
			int row = getSelectedClass();
			gallery.removeClass(row);
			model.fireTableRowsDeleted(row, row);
		}
	}// function actionPerformed

	public int getSelectedClass() {
		return tableClasses.getSelectedRow();
	}

	public ClassInfo getClassInfo(int index) {
		return classes.get(index);
	}

	/** Structure to store class info */
	public class ClassInfo {
		private String comment;
		private ColorIcon icon;
		public int numberOfClasses; //Classes assigned to superclass
		public int numberOfImages; //total images assigned to superclass

		/** Constructor */
		ClassInfo(String name, Color c) {
			comment = name;
			icon = new ColorIcon(c, 16, 16, 3, true);
		}

		public String getComment() {
			return comment;
		}

		public void setComment(String value) {
			comment = value;
		}

		public Color getColor() {
			return icon.getColor();
		}

		public void setColor(Color value) {
			icon.setColor(value);
		}

		public Icon getIcon() {
			return icon;
		}
	}// class ClassInfo

	/** Table model based on the ArrayList of ClassInfo */
	class ClassesTableModel extends AbstractTableModel {
		private static final long serialVersionUID = 1L;

		private String[] columns = { "Comment", "# Classes",
				"# Images" , "Color"};

		public ClassesTableModel() {
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
			return classes.size();
			// return frame.getParticlePicker().getFamilies().size();
		}

		@Override
		public boolean isCellEditable(int row, int column) {
			return (column == 0 || column == 3);
		}

		@Override
		public void setValueAt(Object value, int row, int column) {
			ClassInfo cli = classes.get(row);
			switch (column) {
			case 0:
				cli.setComment((String) value);
				break;
			case 3:
				cli.setColor((Color) value);
				gallery.fireTableDataChanged();
				break;
			default:
				XmippDialog.showError(parent,
						String.format("Column %d is not editable", column));
			}

		}

		@Override
		public Object getValueAt(int row, int column) {
			ClassInfo cli = classes.get(row);
			switch (column) {
			case 0:
				return cli.getComment();
			case 1:
				return cli.numberOfClasses;
			case 2:
				return cli.numberOfImages;
			case 3:
				return cli.getColor();
			default:
				return null;
			}
		}
		
		public void addNewRow(){
			int newPos = classes.size();
			classes.add(new ClassInfo("NewClass", getNextColor()));			
			fireTableRowsInserted(newPos, newPos);
		}
		
		final int amount  = 100;
		final int lowerLimit = 0x101010;
		final int upperLimit = 0xE0E0D0;
		final int colorStep = (upperLimit-lowerLimit)/amount;
		int colorCounter = 0;
		
		Color getNextColor() {
		    return new Color(upperLimit-colorStep*colorCounter++);
		}

	}// class ColumnsTableModel

}// class ColumnsJDialog
