package gui;

import gui.EditFamiliesJDialog.FamiliesTableModel;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.NumberFormat;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JColorChooser;
import javax.swing.JDialog;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;


import model.Family;


public class AddFamilyJDialog extends JDialog implements ActionListener {

	private JButton addbt;
	private JButton cancelbt;
	private JTextField nametf;
	private JButton colorbt;
	private double position = 0.9;
	private Color color;
	private JColorChooser colorChooser;
	EditFamiliesJDialog parent;
	private JSlider sizesl;
	private JFormattedTextField sizetf;

	public AddFamilyJDialog(EditFamiliesJDialog parent, boolean modal) {
		super(parent, modal);
		setResizable(false);
		this.parent = parent;
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Add Family");
		setLayout(new GridBagLayout());
		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(5, 5, 5, 5);
		constraints.anchor = GridBagConstraints.WEST;

		add(new JLabel("Name"),
				WindowUtils.updateConstraints(constraints, 0, 0, 1));
		nametf = new JTextField(20);
		add(nametf, WindowUtils.updateConstraints(constraints, 1, 0, 1));
		add(new JLabel("Color"),
				WindowUtils.updateConstraints(constraints, 0, 1, 1));
		colorbt = new JButton();
		color = Family.getNextColor();
		colorbt.setIcon(new ColorIcon(color));
		colorbt.setBorderPainted(false);
		add(colorbt, WindowUtils.updateConstraints(constraints, 1, 1, 1));

		int size = Family.getDefaultSize();
		JPanel sizepn = new JPanel();
		sizesl = new JSlider(0, 500, size);
		sizepn.add(sizesl);
		sizetf = new JFormattedTextField(NumberFormat.getIntegerInstance());;
		sizetf.setText(Integer.toString(size));
		sizepn.add(sizetf);
		add(new JLabel("Size"),
				WindowUtils.updateConstraints(constraints, 0, 2, 1));
		add(sizepn, WindowUtils.updateConstraints(constraints, 1, 2, 1));
		addbt = new JButton("Add");
		getRootPane().setDefaultButton(addbt);
		cancelbt = new JButton("Cancel");

		add(addbt, WindowUtils.updateConstraints(constraints, 0, 3, 1));
		add(cancelbt, WindowUtils.updateConstraints(constraints, 1, 3, 1));
		setListeners();
		pack();
		WindowUtils.centerScreen(position, this);
		setVisible(true);
	}

	private void setListeners() {

		colorbt.addActionListener(this);

		addbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				String name = nametf.getText();
				if(sizetf.getValue() == null)
				{
					JOptionPane.showMessageDialog(AddFamilyJDialog.this,
							model.Constants.getEmptyFieldMsg("Size"));
					return;
				}
				int size = ((Number)sizetf.getValue()).intValue();
				try {
					Family g = new Family(name, color, size);

					AddFamilyJDialog.this.parent.addGroup(g);
					setVisible(false);
					dispose();
				} catch (IllegalArgumentException ex) {
					JOptionPane.showMessageDialog(AddFamilyJDialog.this,
							ex.getMessage());
				}
				
			}
		});
		cancelbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				setVisible(false);
				dispose();

			}
		});

		sizetf.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				int size = Integer.parseInt(sizetf.getText());
				sizesl.setValue(size);
			}
		});

		sizesl.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {
				int size = sizesl.getValue();
				sizetf.setText(Integer.toString(size));
			}
		});
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		// Set up the dialog that the button brings up.
		colorChooser = new JColorChooser();
		JDialog dialog = JColorChooser.createDialog(colorbt, "Pick a Color",
				true, // modal
				colorChooser, new ActionListener() {

					@Override
					public void actionPerformed(ActionEvent e) {
						AddFamilyJDialog.this.color = AddFamilyJDialog.this.colorChooser
								.getColor();
						AddFamilyJDialog.this.colorbt.setIcon(new ColorIcon(
								color));
					}
				}, // OK button handler
				null); // no CANCEL button handler
		WindowUtils.centerScreen(AddFamilyJDialog.this.position, dialog);
		dialog.setVisible(true);
	}

}
