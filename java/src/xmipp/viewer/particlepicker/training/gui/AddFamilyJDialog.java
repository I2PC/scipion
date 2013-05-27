package xmipp.viewer.particlepicker.training.gui;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.NumberFormat;

import javax.swing.JButton;
import javax.swing.JColorChooser;
import javax.swing.JDialog;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import xmipp.utils.ColorIcon;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.ColorHelper;
import xmipp.viewer.particlepicker.Family;

public class AddFamilyJDialog extends JDialog implements ActionListener {

	private JButton addbt;
	private JButton cancelbt;
	private JTextField nametf;
	private JFormattedTextField templatestf;
	private JButton colorbt;
	private float position = 0.9f;
	private Color color;
	private JColorChooser colorChooser;
	EditFamiliesJDialog parent;
	private JSlider sizesl;
	private JFormattedTextField sizetf;
	private JPanel sizepn;

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
				XmippWindowUtil.getConstraints(constraints, 0, 0, 1));

		nametf = new JTextField(20);

		add(nametf, XmippWindowUtil.getConstraints(constraints, 1, 0, 1));
		add(new JLabel("Color"),
				XmippWindowUtil.getConstraints(constraints, 0, 1, 1));
		colorbt = new JButton();
		color = ColorHelper.getNextColor();
		colorbt.setIcon(new ColorIcon(color));
		colorbt.setBorderPainted(false);
		add(colorbt, XmippWindowUtil.getConstraints(constraints, 1, 1, 1));

		add(new JLabel("Size"),
				XmippWindowUtil.getConstraints(constraints, 0, 2, 1));
		initSizePane();

		add(sizepn, XmippWindowUtil.getConstraints(constraints, 1, 2, 1));

		add(new JLabel("Templates"),
				XmippWindowUtil.getConstraints(constraints, 0, 3, 1));
		templatestf = new JFormattedTextField(NumberFormat.getIntegerInstance());
		templatestf.setValue(1);
		templatestf.setColumns(2);
		add(templatestf, XmippWindowUtil.getConstraints(constraints, 1, 3, 1));

		addbt = XmippWindowUtil.getTextButton("Add", null);
		getRootPane().setDefaultButton(addbt);
		cancelbt = XmippWindowUtil.getTextButton("Cancel", null);

		add(addbt, XmippWindowUtil.getConstraints(constraints, 0, 4, 1));
		add(cancelbt, XmippWindowUtil.getConstraints(constraints, 1, 4, 1));
		setListeners();
		pack();
		XmippWindowUtil.setLocation(position, 0.5f, this);
		setVisible(true);
	}

	private void setListeners() {

		colorbt.addActionListener(this);

		addbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				String name = nametf.getText();
				if (sizetf.getValue() == null)
					sizetf.setValue(sizesl.getValue());
				try {
					if (templatestf.getValue() == null || templatestf.getText().equals(""))
						throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("Templates"));
					int size = ((Number) sizetf.getValue()).intValue();
					int templatesNumber = ((Number) templatestf.getValue()).intValue();
					if(templatesNumber < 1)
						throw new IllegalArgumentException(XmippMessage.getIllegalValueMsg("Templates", templatesNumber));
					Family g = new Family(name, color, size, parent.getFrame().getParticlePicker(), templatesNumber);

					AddFamilyJDialog.this.parent.addFamily(g);
					setVisible(false);
					dispose();
				} catch (IllegalArgumentException ex) {
					ex.printStackTrace();
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

	}

	private void initSizePane() {
		int size = Family.getDefaultSize();
		sizepn = new JPanel();
		sizesl = new JSlider(0, 500, size);
		sizesl.setMajorTickSpacing(250);
		sizesl.setPaintTicks(true);
		sizesl.setPaintLabels(true);
		sizepn.add(sizesl);
		sizetf = new JFormattedTextField(NumberFormat.getIntegerInstance());
		
		sizetf.setValue(size);
		sizepn.add(sizetf);
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
				sizetf.setValue(size);
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
		XmippWindowUtil.setLocation(AddFamilyJDialog.this.position, 0.5f,
				dialog);
		dialog.setVisible(true);
	}

}
