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
import xmipp.utils.ColorRenderer;
import xmipp.utils.XmippWindowUtil;
import xmipp.utils.XmippDialog;
import xmipp.viewer.models.ColumnInfo;

public class ClassesJDialog extends XmippDialog {
	private static final long serialVersionUID = 1L;
	private JTable tableClasses;
	private JButton btnUp;
	private JButton btnDown;
	private ClassesTableModel model;
	// This will be used for check for results from the dialog
	private ArrayList<ClassInfo> classes;
	boolean fireEvent = true;

	public ClassesJDialog(JFrameGallery parent) {
		super(parent, "Columns", true);
		disposeOnClose = false;
		initComponents();
	}// constructor ColumnsJDialog

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
				.createTitledBorder("Classes"));
		groupstbpn.add(sp);
		sp.setOpaque(true);
		model = new ClassesTableModel();
		tableClasses = new JTable(model);
		tableClasses
				.setPreferredScrollableViewportSize(new Dimension(350, 200));
		sp.setViewportView(tableClasses);
		tableClasses.setDefaultRenderer(Color.class, new ColorRenderer());
		tableClasses.setDefaultEditor(Color.class, new ColorEditor());
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
		tableClasses.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		tableClasses.getSelectionModel().addListSelectionListener(
				new ListSelectionListener() {
					@Override
					public void valueChanged(ListSelectionEvent e) {
						enableUpDown(true);
					}
				});
	}// function initComponents

	private void enableUpDown(boolean value) {
		btnUp.setEnabled(value);
		btnDown.setEnabled(value);
	}// function enableUpDown

	@Override
	public void handleActionPerformed(ActionEvent evt){		
		JButton btn = (JButton) evt.getSource();

	}// function actionPerformed

	public int getSelectedClass(){
		return tableClasses.getSelectedRow();
	}
	
	public ClassInfo getClassInfo(int index){
		return classes.get(index);
	}
	
	/** Structure to store class info */
	public class ClassInfo {
		private String comment;
		private Color color;
		
		/** Constructor */
		ClassInfo(String name, Color c){
			comment = name;
			color = c;
		}
		
		public String getComment(){
			return comment;
		}
		
		public void setComment(String value){
			comment = value;
		}
		
		public Color getColor(){
			return color;
		}
		
		public void setColor(Color value){
			color = value;
		}
	}//class ClassInfo
	
	
	/** Table model based on the ArrayList of ClassInfo */
	class ClassesTableModel extends AbstractTableModel {
		private static final long serialVersionUID = 1L;

		private String[] columns = { "Comment", "Color", "# Classes", "# Images" };

		public ClassesTableModel() {
			classes = new ArrayList<ClassInfo>();
			classes.add(new ClassInfo("kk", Color.red));
			classes.add(new ClassInfo("otra clase", Color.blue));
			classes.add(new ClassInfo("kk", Color.green));
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
			return (column <= 1);
		}

		@Override
		public void setValueAt(Object value, int row, int column) {
			ClassInfo cli = classes.get(row);
			switch (column){
			case 0:
				cli.setComment((String)value);
				break;
			case 1:
				cli.setColor((Color)value);
				break;
			default:
				XmippDialog.showError(parent, 
						String.format("Column %d is not editable", column));
			}

		}

		@Override
		public Object getValueAt(int row, int column) {
			ClassInfo cli = classes.get(row);
			switch (column){
			case 0:
				return cli.getComment();
			case 1:
				return cli.getColor();
			default:
				return 0;
			}
		}

	}// class ColumnsTableModel

}// class ColumnsJDialog
