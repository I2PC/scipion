package xmipp.particlepicker.training.gui;


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

import xmipp.particlepicker.Family;
import xmipp.utils.WindowUtils;




public class AddFamilyJDialog extends JDialog implements ActionListener {

	private JButton addbt;
	private JButton cancelbt;
	private JTextField nametf;
	private JButton colorbt;
	private float position = 0.9f;
	private Color color;
	private JColorChooser colorChooser;
	EditFamiliesJDialog parent;
	private JSlider sizesl;
	private JFormattedTextField sizetf;
	private JPanel sizepn;
//	private JPanel thresholdpn;
//	private JSlider thresholdsl;
//	private JFormattedTextField thresholdtf;
	private int range;

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
				WindowUtils.getConstraints(constraints, 0, 0, 1));
		nametf = new JTextField(20);
		add(nametf, WindowUtils.getConstraints(constraints, 1, 0, 1));
		add(new JLabel("Color"),
				WindowUtils.getConstraints(constraints, 0, 1, 1));
		colorbt = new JButton();
		color = Family.getNextColor();
		colorbt.setIcon(new ColorIcon(color));
		colorbt.setBorderPainted(false);
		add(colorbt, WindowUtils.getConstraints(constraints, 1, 1, 1));

		
		add(new JLabel("Size"),
				WindowUtils.getConstraints(constraints, 0, 2, 1));
		initSizePane();
		
		add(sizepn, WindowUtils.getConstraints(constraints, 1, 2, 1));
//		add(new JLabel("Threshold"),
//				WindowUtils.updateConstraints(constraints, 0, 3, 1));
//		initThresholdPane();
//		add(thresholdpn, WindowUtils.updateConstraints(constraints, 1, 3, 1));
		addbt = new JButton("Add");
		getRootPane().setDefaultButton(addbt);
		cancelbt = new JButton("Cancel");

		add(addbt, WindowUtils.getConstraints(constraints, 0, 3, 1));
		add(cancelbt, WindowUtils.getConstraints(constraints, 1, 3, 1));
		setListeners();
		pack();
		WindowUtils.setLocation(position, 0.5f, this);
		setVisible(true);
	}


	private void setListeners() {

		colorbt.addActionListener(this);

		addbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				String name = nametf.getText();
				if(sizetf.getValue() == null)
					sizetf.setValue(sizesl.getValue());
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

		
	}
	
	private void initSizePane()
	{
		int size = Family.getDefaultSize();
		sizepn = new JPanel();
		sizesl = new JSlider(0, 500, size);
		sizesl.setMajorTickSpacing(250);
		sizesl.setPaintTicks(true);
		sizesl.setPaintLabels(true);
		sizepn.add(sizesl);
		sizetf = new JFormattedTextField(NumberFormat.getIntegerInstance());;
		sizetf.setText(Integer.toString(size));
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
				sizetf.setText(Integer.toString(size));
			}
		});
	}
	
//	private void initThresholdPane() {
//		int threshold = 0;
//		range = 10;
//		thresholdpn = new JPanel();
//		thresholdsl = new JSlider(0, range);
//		Hashtable<Integer, JComponent> labelTable = new Hashtable<Integer, JComponent>();
//		labelTable.put( new Integer( 0 ), new JLabel("0.0") );
//		labelTable.put( new Integer( range/4 ), new JLabel("0.25") );
//		labelTable.put( new Integer( range/2 ), new JLabel("0.5") );
//		labelTable.put( new Integer( 3*range/4 ), new JLabel("0.75") );
//		labelTable.put( new Integer( range ), new JLabel("1.0") );
//		thresholdpn.add(thresholdsl);
//		thresholdtf = new JFormattedTextField(NumberFormat.getNumberInstance());;
//		thresholdtf.setColumns(3);
//		thresholdtf.setText(Integer.toString(threshold));
//		thresholdpn.add(thresholdtf);
//		thresholdtf.addActionListener(new ActionListener() {
//
//			@Override
//			public void actionPerformed(ActionEvent e) {
//				double threshold = Double.parseDouble(thresholdtf.getText() ) * range;
//				int range = AddFamilyJDialog.this.range;
//				if(Math.abs(threshold) <= range)
//					thresholdsl.setValue((int)threshold );
//			}
//		});
//
//		thresholdsl.addChangeListener(new ChangeListener() {
//
//			@Override
//			public void stateChanged(ChangeEvent e) {
//				int range = AddFamilyJDialog.this.range;
//				double threshold = (double)thresholdsl.getValue()/range;
//				thresholdtf.setText(String.format("%.2f", threshold));
//			}
//		});
//		
//	}


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
		WindowUtils.setLocation(AddFamilyJDialog.this.position, 0.5f, dialog);
		dialog.setVisible(true);
	}

}
