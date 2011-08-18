package ppicker;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;

public class EditFamiliesJDialog extends JDialog {

	private XmippParticlePickerJFrame parent;
	private JTable groupstb;
	private JButton addbt;
	private JButton deletebt;
	private JButton okbt;
	private FamiliesTableModel model;

	public EditFamiliesJDialog(XmippParticlePickerJFrame parent, boolean modal) {
		super(parent, modal);
		this.parent = parent;
		initComponents();
	}

	private void initComponents() {
		setResizable(false);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Families");
		setLayout(new GridBagLayout());
		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(10, 10, 10, 10);
		constraints.anchor = GridBagConstraints.WEST;

		JPanel groupstbpn = new JPanel();
		JScrollPane sp = new JScrollPane();
		groupstbpn.setBorder(BorderFactory.createTitledBorder("Families"));
		groupstbpn.add(sp);
		sp.setOpaque(true);
		double position = 0.9;
		model = new FamiliesTableModel(parent);
		groupstb = new JTable(model);
		groupstb.setPreferredScrollableViewportSize(new Dimension(200, 200));
		groupstb.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		groupstb.setDefaultRenderer(Color.class, new ColorRenderer());
		groupstb.setDefaultEditor(Color.class, new ColorEditor(position));
		sp.setViewportView(groupstb);
		add(groupstbpn, WindowUtils.updateConstraints(constraints, 0, 0, 3));
		addbt = new JButton("Add");
		deletebt = new JButton("Delete");
		deletebt.setEnabled(false);
		okbt = new JButton("Ok");

		add(addbt, WindowUtils.updateConstraints(constraints, 0, 1, 1));
		add(deletebt, WindowUtils.updateConstraints(constraints, 1, 1, 1));
		add(okbt, WindowUtils.updateConstraints(constraints, 2, 1, 1));
		getRootPane().setDefaultButton(addbt);
		setListeners();
		pack();
		WindowUtils.centerScreen(position, this);
		setVisible(true);
	}

	private void setListeners() {
		
		addbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				new AddFamilyJDialog(EditFamiliesJDialog.this, true);
			}
		});

		groupstb.getSelectionModel().addListSelectionListener(
				new ListSelectionListener() {

					@Override
					public void valueChanged(ListSelectionEvent e) {
						EditFamiliesJDialog.this.deletebt.setEnabled(true);

					}
				});
		okbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				setVisible(false);
				dispose();

			}
		});
		deletebt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				List<Family> groups = EditFamiliesJDialog.this.parent
						.getPPData().getFamilies();
				if (groups.size() == 1) {
					JOptionPane.showMessageDialog(EditFamiliesJDialog.this,
							Constants.delete_last_group_msg);
					return;
				}
				int index = groupstb.getSelectedRow();
				groups.remove(index);
				model.fireTableStructureChanged();
				EditFamiliesJDialog.this.deletebt.setEnabled(false);
				parent.updateFamilies();
			}
		});
	}

	class FamiliesTableModel extends AbstractTableModel {

		private String[] columns = new String[] { "Name", "Color", "Radius" };
		private XmippParticlePickerJFrame frame;

		public FamiliesTableModel(XmippParticlePickerJFrame frame) {

			this.frame = frame;
		}
		
		@Override
		public Class getColumnClass(int column)
		{
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
			return frame.getPPData().getFamilies().size();
		}

		@Override
		public boolean isCellEditable(int row, int column) {
			return true;
		}

		@Override
		public void setValueAt(Object value, int row, int column) {
			Family g = frame.getPPData().getFamilies().get(row);
			if (column == 0) {
				String name = (String) value;
				if(name.equals(g.getName()))
						return;
				if (name.equals("")) {
					JOptionPane.showMessageDialog(EditFamiliesJDialog.this,
							Constants.empty_group_name_msg);
					return;
				}
				else if(parent.getPPData().existsGroupName(name))
					JOptionPane.showMessageDialog(EditFamiliesJDialog.this, Constants.getAlreadyExistsGroupNameMsg(name));
				g.setName(name);
			} else if (column == 1)
				g.setColor((Color) value);
			else if(column == 2)
			{
				int radius = (Integer)value;
				g.setRadius(radius);
			}
			frame.updateFamilies();
		}

		@Override
		public Object getValueAt(int row, int column) {
			Family f = frame.getPPData().getFamilies().get(row);
			if (column == 0)
				return f.getName();
			if (column == 1)
				return f.getColor();
			if(column == 2)
				return f.getRadius();
			return null;
		}

	}

	public void addGroup(Family g) {
		parent.addGroup(g);
		model.fireTableStructureChanged();
	}

}
