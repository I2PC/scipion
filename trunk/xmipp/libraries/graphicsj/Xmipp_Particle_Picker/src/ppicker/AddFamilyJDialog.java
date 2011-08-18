package ppicker;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JColorChooser;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import ppicker.EditFamiliesJDialog.FamiliesTableModel;

public class AddFamilyJDialog extends JDialog implements ActionListener {

	private JButton addbt;
	private JButton cancelbt;
	private JTextField nametf;
	private JButton colorbt;
	private double position = 0.9;
	private Color color;
	private JColorChooser colorChooser;
	EditFamiliesJDialog parent;
	private JSlider radiussl;
	private JTextField radiustf;

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
		Color[] colors = Family.getSampleColors();
		int index = (int)(Math.random() * colors.length);
		color = colors[index];
		colorbt.setIcon(new ColorIcon(color));
		colorbt.setBorderPainted(false);
		add(colorbt, WindowUtils.updateConstraints(constraints, 1, 1, 1));

		int radius = Family.getDefaultRadius();
		JPanel radiuspn = new JPanel();
		radiussl = new JSlider(0, 500, radius);
		radiuspn.add(radiussl);
		radiustf = new JTextField(3);
		radiustf.setText(Integer.toString(radius));
		radiuspn.add(radiustf);
		add(new JLabel("Radius"), WindowUtils.updateConstraints(constraints, 0, 2, 1));
		add(radiuspn, WindowUtils.updateConstraints(constraints, 1, 2, 1));
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
				if (name.equals("")) {
					JOptionPane.showMessageDialog(AddFamilyJDialog.this,
							Constants.empty_group_name_msg);
					return;
				}
				int radius = (Integer)radiussl.getValue();
				Family g = new Family(name, color, radius);
				try
				{
					AddFamilyJDialog.this.parent.addGroup(g);
				}
				catch(IllegalArgumentException ex)
				{
					JOptionPane.showMessageDialog(AddFamilyJDialog.this, ex.getMessage());
				}
				setVisible(false);
				dispose();
			}
		});
		cancelbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				setVisible(false);
				dispose();

			}
		});
		
		radiustf.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				int radius = Integer.parseInt(radiustf.getText());
				radiussl.setValue(radius);
			}
		});

		radiussl.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {
				int radius = radiussl.getValue();
				radiustf.setText(Integer.toString(radius));
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
						AddFamilyJDialog.this.colorbt.setIcon(new ColorIcon(color));
					}
				}, // OK button handler
				null); // no CANCEL button handler
		WindowUtils.centerScreen(AddFamilyJDialog.this.position, dialog);
		dialog.setVisible(true);
	}

}
