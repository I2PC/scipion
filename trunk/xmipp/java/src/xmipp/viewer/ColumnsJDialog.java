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
import javax.swing.table.AbstractTableModel;

import xmipp.jni.MetaData;
import xmipp.utils.WindowUtil;



public class ColumnsJDialog extends JDialog {

	private JFrameGallery parent;
	private JTable groupstb;
	private JButton cancelbt;
	private JButton okbt;
	private ColumnsTableModel model;
	//This will be used for check for results from the dialog
	private ArrayList<ColumnInfo> rows;

	public ColumnsJDialog(JFrameGallery parent, boolean modal) {
		super(parent, modal);
		this.parent = parent;
		initComponents();
	}

	public ArrayList<ColumnInfo> getColumnsResult(){
		return rows;
	}
	
	private void initComponents() {
		setResizable(false);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Columns");
		setLayout(new GridBagLayout());
		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(10, 10, 10, 10);
		constraints.anchor = GridBagConstraints.WEST;

		JPanel groupstbpn = new JPanel();
		JScrollPane sp = new JScrollPane();
		groupstbpn.setBorder(BorderFactory.createTitledBorder("Column properties"));
		groupstbpn.add(sp);
		sp.setOpaque(true);
		model = new ColumnsTableModel(parent.getData().labels);
		groupstb = new JTable(model);
		groupstb.setPreferredScrollableViewportSize(new Dimension(350, 200));
		groupstb.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		sp.setViewportView(groupstb);		
		add(groupstbpn, WindowUtil.getConstraints(constraints, 0, 0, 3));		
		cancelbt = new JButton("Cancel");
		okbt = new JButton("Ok");
		add(okbt, WindowUtil.getConstraints(constraints, 1, 1, 1));
		add(cancelbt, WindowUtil.getConstraints(constraints, 2, 1, 1));
		getRootPane().setDefaultButton(okbt);
		
		addListeners();
		pack();
		WindowUtil.centerWindows(this, parent);
		setVisible(true);
	}

	private void addListeners() {

		okbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				setVisible(false);
				dispose();

			}
		});
		cancelbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				rows = null; 
				setVisible(false);
				dispose();
			}
		});
	}

	class ColumnsTableModel extends AbstractTableModel {
		private static final long serialVersionUID = 1L;

		private String[] columns = {"Label", "Visible", "Render"};
		
		
		public ColumnsTableModel(int[] labels){
			rows = new ArrayList<ColumnInfo>(labels.length);
			for (int i = 0; i < labels.length; ++i)
				rows.add(new ColumnInfo(labels[i]));
		}
		
		public ColumnsTableModel(ArrayList<ColumnInfo> labelsInfo){
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
			//return frame.getParticlePicker().getFamilies().size();
		}

		@Override
		public boolean isCellEditable(int row, int column) {
			try {
				if (column == 0 || 
				   (column == 2 && !MetaData.isImage(rows.get(row).getLabel())))
					return false;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return true;
		}

		@Override
		public void setValueAt(Object value, int row, int column) {
			ColumnInfo col = rows.get(row);
			switch (column){
			case 1:
				col.visible = (Boolean)value;
				break;
			case 2:
				col.render = (Boolean)value;
				break;
			}
//			try {
//				if(value == null)
//					throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg(columns[column]));
//				Family f = frame.getParticlePicker().getFamilies().get(row);
//				if (column == 0) {
//					String name = (String) value;
//					if (name.equals(f.getName()))
//						return;
//					else if (parent.getParticlePicker().existsFamilyName(name))
//						JOptionPane.showMessageDialog(ColumnsJDialog.this,
//								XmippMessage.getAlreadyExistsGroupNameMsg(name));
//					f.setName(name);
//					frame.updateFamilyComboBox();
//				} else if (column == 1)
//				{
//					f.setColor((Color) value);
//					frame.updateFamilyColor();
//				}
//				else if (column == 2) {
//					int size = (Integer) value;
//					f.setSize(size);
//					frame.switchSize(size);
//				}
//
//			} catch (IllegalArgumentException e) {
//				JOptionPane.showMessageDialog(ColumnsJDialog.this,
//						e.getMessage());
//			}
		}

		@Override
		public Object getValueAt(int row, int column) {
			ColumnInfo col = rows.get(row);
			switch (column){
			case 0:
				return col.getLabelName();
			case 1:
				return col.visible;
			case 2:
				return col.render;
			}
			return null;
		}

	}

}
