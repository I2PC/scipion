package xmipp.viewer.particlepicker.training.gui;

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

import xmipp.utils.ColorEditor;
import xmipp.utils.ColorRenderer;
import xmipp.utils.XmippWindowUtil;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.Family;

public class EditFamiliesJDialog extends JDialog {

	private TrainingPickerJFrame parent;
	private JTable familiestb;
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
		float position = 0.9f;
		model = new FamiliesTableModel(parent);
		familiestb = new JTable(model);
		familiestb.setPreferredScrollableViewportSize(new Dimension(300, 200));
		familiestb
				.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		familiestb.setDefaultRenderer(Color.class, new ColorRenderer());
		familiestb.setDefaultEditor(Color.class, new ColorEditor(position));
		sp.setViewportView(familiestb);
		add(groupstbpn, XmippWindowUtil.getConstraints(constraints, 0, 0, 3));
		addbt = XmippWindowUtil.getTextButton("Add", null);
		deletebt = XmippWindowUtil.getTextButton("Delete", null);
		deletebt.setEnabled(false);
		okbt = XmippWindowUtil.getTextButton("Ok", null);

		add(addbt, XmippWindowUtil.getConstraints(constraints, 0, 1, 1));
		add(deletebt, XmippWindowUtil.getConstraints(constraints, 1, 1, 1));
		add(okbt, XmippWindowUtil.getConstraints(constraints, 2, 1, 1));
		getRootPane().setDefaultButton(addbt);
		setListeners();
		pack();
		XmippWindowUtil.setLocation(position, 0.5f, this);
		setVisible(true);
	}

	private void setListeners() {

		addbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				new AddFamilyJDialog(EditFamiliesJDialog.this, true);
			}
		});

		familiestb.getSelectionModel().addListSelectionListener(
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
				List<Family> families = parent.getParticlePicker()
						.getFamilies();
				try {

					int index = familiestb.getSelectedRow();
					Family f = families.get(index);
					if (f.equals(parent.getFamily()))
						throw new IllegalArgumentException(
								"Can not remove active family");
					parent.removeFamily(f);
					model.fireTableStructureChanged();
					EditFamiliesJDialog.this.deletebt.setEnabled(false);
				} catch (Exception ex) {
					JOptionPane.showMessageDialog(EditFamiliesJDialog.this,
							ex.getMessage());

				}
			}
		});
	}

	class FamiliesTableModel extends AbstractTableModel {

		private String[] columns = new String[] { "Name", "Color", "Size", "Templates" };
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
				if (value == null)
					throw new IllegalArgumentException(
							XmippMessage.getEmptyFieldMsg(columns[column]));
				Family f = frame.getParticlePicker().getFamilies().get(row);
				if (column == 0) {
					String name = (String) value;
					if (name.equals(f.getName()))
						return;
					else if (parent.getParticlePicker().existsFamilyName(name)) {
						JOptionPane
								.showMessageDialog(
										EditFamiliesJDialog.this,
										XmippMessage
												.getAlreadyExistsGroupNameMsg(name));
						return;
					}

					f.setName(name);
					frame.updateFamilyComboBox();
				} else if (column == 1) {
					f.setColor((Color) value);
					frame.updateFamilyColor();
				} else if (column == 2) {
					int size = (Integer) value;
					frame.updateSize(size);
				
				} else if (column == 3) {
					
					int templates = (Integer)value;
					if(templates  < 1)
						throw new IllegalArgumentException(XmippMessage.getIllegalValueMsg("Templates", templates));
					frame.setTemplatesNumber(f, templates);
				}
				frame.getParticlePicker().saveFamilies();
				

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
			if (column == 3)
				return f.getTemplatesNumber();
			return null;
		}

	}

	public void addFamily(Family g) {

		parent.addFamily(g);
		model.fireTableStructureChanged();
	}

	public TrainingPickerJFrame getFrame()
	{
		return parent;
	}

}
