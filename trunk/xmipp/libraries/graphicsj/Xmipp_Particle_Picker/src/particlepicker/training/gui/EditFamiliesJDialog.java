package particlepicker.training.gui;

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

import particlepicker.Constants;
import particlepicker.Family;
import particlepicker.WindowUtils;



public class EditFamiliesJDialog extends JDialog {

	private TrainingPickerJFrame parent;
	private JTable groupstb;
	private JButton addbt;
	private JButton deletebt;
	private JButton okbt;
	private FamiliesTableModel model;

	public EditFamiliesJDialog(TrainingPickerJFrame parent, boolean modal) {
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
		add(groupstbpn, WindowUtils.getConstraints(constraints, 0, 0, 3));
		addbt = new JButton("Add");
		deletebt = new JButton("Delete");
		deletebt.setEnabled(false);
		okbt = new JButton("Ok");

		add(addbt, WindowUtils.getConstraints(constraints, 0, 1, 1));
		add(deletebt, WindowUtils.getConstraints(constraints, 1, 1, 1));
		add(okbt, WindowUtils.getConstraints(constraints, 2, 1, 1));
		getRootPane().setDefaultButton(addbt);
		setListeners();
		pack();
		WindowUtils.centerScreen(position, 0.5, this);
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
				List<Family> families = EditFamiliesJDialog.this.parent
						.getParticlePicker().getFamilies();
				try
				{
				int index = groupstb.getSelectedRow();
				parent.removeFamily(families.get(index));
				model.fireTableStructureChanged();
				EditFamiliesJDialog.this.deletebt.setEnabled(false);
				}
				catch(Exception ex)
				{
					JOptionPane.showMessageDialog(EditFamiliesJDialog.this,
							ex.getMessage());
				
				}
			}
		});
	}

	class FamiliesTableModel extends AbstractTableModel {

		private String[] columns = new String[] { "Name", "Color", "Size" };
		private TrainingPickerJFrame frame;

		public FamiliesTableModel(TrainingPickerJFrame frame) {

			this.frame = frame;
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
			return frame.getParticlePicker().getFamilies().size();
		}

		@Override
		public boolean isCellEditable(int row, int column) {
			return true;
		}

		@Override
		public void setValueAt(Object value, int row, int column) {
			try {
				if(value == null)
					throw new IllegalArgumentException(Constants.getEmptyFieldMsg(columns[column]));
				Family f = frame.getParticlePicker().getFamilies().get(row);
				if (column == 0) {
					String name = (String) value;
					if (name.equals(f.getName()))
						return;
					else if (parent.getParticlePicker().existsFamilyName(name))
						JOptionPane.showMessageDialog(EditFamiliesJDialog.this,
								Constants.getAlreadyExistsGroupNameMsg(name));
					f.setName(name);
					frame.updateFamilyComboBox();
				} else if (column == 1)
				{
					f.setColor((Color) value);
					frame.updateFamilyColor();
				}
				else if (column == 2) {
					int size = (Integer) value;
					f.setSize(size);
					frame.switchSize(size);
				}

			} catch (IllegalArgumentException e) {
				JOptionPane.showMessageDialog(EditFamiliesJDialog.this,
						e.getMessage());
			}
		}

		@Override
		public Object getValueAt(int row, int column) {
			Family f = frame.getParticlePicker().getFamilies().get(row);
			if (column == 0)
				return f.getName();
			if (column == 1)
				return f.getColor();
			if (column == 2)
				return f.getSize();
			return null;
		}

	}

	public void addGroup(Family g) {
		parent.addFamily(g);
		model.fireTableStructureChanged();
	}

}
