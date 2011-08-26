package gui;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.ImageWindow;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JColorChooser;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import model.State;
import model.Constants;
import model.Family;
import model.Micrograph;
import model.MicrographFamilyData;
import model.Particle;
import model.ParticlePicker;
import model.Step;
import model.XmippJ;
import browser.windows.ImagesWindowFactory;

enum Tool {
	IMAGEJ, PICKER
}

enum Shape {
	Circle, Rectangle, Center
}


public class ParticlePickerJFrame extends JFrame implements ActionListener {

	private JSlider sizesl;
	private JCheckBox circlechb;
	private JCheckBox rectanglechb;
	private ParticlePickerCanvas canvas;
	private JFormattedTextField sizetf;
	private Tool tool = Tool.PICKER;
	private JCheckBox centerchb;
	private JMenuBar mb;
	private JComboBox familiescb;
	private ParticlePicker ppicker;
	private Color color;
	private JPanel familypn;
	private JPanel symbolpn;
	private String activemacro;
	private JPanel micrographpn;
	private JTable micrographstb;
	private ImageWindow iw;
	private boolean changed;
	private JMenuItem savemi;
	private MicrographsTableModel micrographsmd;
	Micrograph micrograph;
	private JButton nextbt;
	private JButton colorbt;
	private double position;
	private JLabel iconlb;
	private JLabel steplb;
	private JButton actionsbt;
	private Family family;
	private JMenuItem editfamiliesmn;

	// private JCheckBox onlylastchb;

	public boolean isShapeSelected(Shape s) {
		switch (s) {
		case Rectangle:
			return rectanglechb.isSelected();
		case Circle:
			return circlechb.isSelected();
		case Center:
			return centerchb.isSelected();
			// case OnlyLast:
			// return onlylastchb.isSelected();
		}
		return false;
	}

	public Tool getTool() {
		return tool;
	}

	public ParticlePicker getParticlePicker() {
		return ppicker;
	}

	public Family getFamily() {
		return family;
	}

	public ParticlePickerJFrame() {

		ppicker = ParticlePicker.getInstance();

		initComponents();
		initializeCanvas();
	}

	public Micrograph getMicrograph() {
		return micrograph;
	}

	private void initComponents() {
		// try {
		// UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		// } catch (Exception e) {
		// // TODO Auto-generated catch block
		// e.printStackTrace();
		// }

		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent winEvt) {
				if (changed) {
					int result = JOptionPane.showConfirmDialog(
							ParticlePickerJFrame.this,
							"Save changes before closing?", "Message",
							JOptionPane.YES_NO_OPTION);
					if (result == JOptionPane.OK_OPTION)
						ParticlePickerJFrame.this.saveChanges();
				}
				System.exit(0);
			}

		});
		setResizable(false);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Particle Picker");
		initPPMenuBar();
		setJMenuBar(mb);

		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(0, 5, 0, 5);
		constraints.anchor = GridBagConstraints.WEST;
		setLayout(new GridBagLayout());

		initFamilyPane();
		add(familypn, WindowUtils.updateConstraints(constraints, 0, 1, 3));

		initSymbolPane();
		add(symbolpn, WindowUtils.updateConstraints(constraints, 0, 2, 3));

		initMicrographsPane();
		add(micrographpn, WindowUtils.updateConstraints(constraints, 0, 3, 3));

		pack();
		position = 0.9;
		WindowUtils.centerScreen(position, this);
		setVisible(true);
	}

	public void initPPMenuBar() {
		mb = new JMenuBar();

		// Setting menus
		JMenu filemn = new JMenu("File");
		JMenu filtersmn = new JMenu("Filters");
		JMenu windowmn = new JMenu("Window");
		JMenu helpmn = new JMenu("Help");
		mb.add(filemn);
		mb.add(filtersmn);
		mb.add(windowmn);
		mb.add(helpmn);

		// Setting menu items
		savemi = new JMenuItem("Save");
		savemi.setEnabled(false);
		filemn.add(savemi);
		JMenuItem stackmi = new JMenuItem("Generate Stack...");
		filemn.add(stackmi);

		JMenuItem bcmi = new JMenuItem("Brightness/Contrast...");
		filtersmn.add(bcmi);
		bcmi.addActionListener(this);

		JMenuItem fftbpf = new JMenuItem("Bandpass Filter...");
		filtersmn.add(fftbpf);
		fftbpf.addActionListener(this);
		JMenuItem admi = new JMenuItem("Anisotropic Diffusion...");
		filtersmn.add(admi);
		admi.addActionListener(this);
		JMenuItem msmi = new JMenuItem("Mean Shift");
		filtersmn.add(msmi);
		msmi.addActionListener(this);
		JMenuItem sbmi = new JMenuItem("Substract Background...");
		filtersmn.add(sbmi);
		sbmi.addActionListener(this);
		JMenuItem gbmi = new JMenuItem("Gaussian Blur...");
		filtersmn.add(gbmi);
		gbmi.addActionListener(this);

		JMenuItem particlesmn = new JMenuItem("Particles");
		windowmn.add(particlesmn);
		JMenuItem ijmi = new JMenuItem("ImageJ");
		windowmn.add(ijmi);
		editfamiliesmn = new JMenuItem("Edit Families");
		windowmn.add(editfamiliesmn);

		JMenuItem hcontentsmi = new JMenuItem("Help Contents...");
		helpmn.add(hcontentsmi);

		// Setting menu item listeners

		ijmi.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				if (IJ.getInstance() == null) {
					new ImageJ(ImageJ.EMBEDDED);
					IJ.getInstance();
				}
				IJ.getInstance().setVisible(true);
			}
		});
		savemi.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				saveChanges();
				JOptionPane.showMessageDialog(ParticlePickerJFrame.this,
						"Data saved successfully");
				((JMenuItem) e.getSource()).setEnabled(false);
			}
		});
		stackmi.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {

			}
		});

		particlesmn.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				List<ImagePlus> imgs = new ArrayList<ImagePlus>();
				for (Particle p : getMicrograph().getFamilyData(family)
						.getManualParticles())
					imgs.add(p.getImage());
				String filename = XmippJ.saveTempImageStack(imgs);
				ImagesWindowFactory.openFileAsImage(filename);
				// new
				// MicrographParticlesJDialog(XmippParticlePickerJFrame.this,
				// XmippParticlePickerJFrame.this.micrograph);
			}
		});

		editfamiliesmn.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				new EditFamiliesJDialog(ParticlePickerJFrame.this, true);

			}
		});

	}

	private void initFamilyPane() {
		familypn = new JPanel();
		familypn.setLayout(new GridLayout(3, 1));
		familypn.setBorder(BorderFactory.createTitledBorder("Family"));

		JPanel fieldspn = new JPanel(new FlowLayout(FlowLayout.LEFT));

		// Setting combo
		fieldspn.add(new JLabel("Family:"));
		familiescb = new JComboBox(ppicker.getFamilies().toArray());
		
		family = (Family) familiescb.getSelectedItem();
		
		fieldspn.add(familiescb);

		// Setting color
		color = family.getColor();
		fieldspn.add(new JLabel("Color:"));
		colorbt = new JButton();
		colorbt.setIcon(new ColorIcon(color));
		colorbt.setBorderPainted(false);
		fieldspn.add(colorbt);

		// Setting slider
		int size = family.getSize();
		fieldspn.add(new JLabel("Size:"));
		sizesl = new JSlider(0, 500, size);
		int height = (int) sizesl.getPreferredSize().getHeight();
		sizesl.setPreferredSize(new Dimension(100, height));

		fieldspn.add(sizesl);
		sizetf = new JFormattedTextField(NumberFormat.getNumberInstance());
		sizetf.setText(Integer.toString(size));
		sizetf.setColumns(3);
		fieldspn.add(sizetf);

		familypn.add(fieldspn, 0);
		JPanel steppn = new JPanel(new FlowLayout(FlowLayout.LEFT));
		steppn.add(new JLabel("Step:"));
		Step step = family.getStep();
		steplb = new JLabel();
		steppn.add(steplb);
		nextbt = new JButton();
		steppn.add(nextbt);
		micrograph = ppicker.getMicrographs().get(0);
		
		JPanel buttonspn = new JPanel(new FlowLayout(FlowLayout.RIGHT));
		actionsbt = new JButton(getFamilyData().getAction());
		actionsbt.setVisible(getFamilyData().isActionAvailable());
		setStep(step, false);
		actionsbt.setText(getFamilyData().getState().toString());
		actionsbt.setVisible(getFamilyData().isActionAvailable());
		buttonspn.add(actionsbt);
		familypn.add(buttonspn);


		colorbt.addActionListener(new ColorActionListener());
		nextbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				try {
					family.validateNextStep();
				} catch (Exception ex) {
					JOptionPane.showMessageDialog(ParticlePickerJFrame.this,
							ex.getMessage());
					return;
				}
				if (family.getStep() == Step.Manual) {
					int result = JOptionPane.showConfirmDialog(
									ParticlePickerJFrame.this,
									"Data provided from manual picking will be used to train software and start Supervised mode."
											+ "\nManual mode will be disabled. Are you sure you want to continue?",
									"Message", JOptionPane.YES_NO_OPTION);

					if (result == JOptionPane.NO_OPTION)
						return;
					setStep(Step.Supervised, true);
					saveChanges();
					ppicker.train(family);
					int next = ppicker.getNextFreeMicrograph(family);
					if(next != -1 )
						micrographstb.setRowSelectionInterval(next, next);
					actionsbt.setText(getFamilyData().getState().toString());
					actionsbt.setVisible(getFamilyData().isActionAvailable());
					
				} else if (family.getStep() == Step.Supervised) {
					int result = JOptionPane.showConfirmDialog(
							ParticlePickerJFrame.this,
							"Model generated during Supervised mode will be dismissed."
									+ "\nAre you sure you want to continue?",
							"Message", JOptionPane.YES_NO_OPTION);

					if (result == JOptionPane.NO_OPTION)
						return;
					setStep(Step.Manual, true);
				}
			}
		});
		familypn.add(steppn, 1);
		sizetf.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				int size = ((Number) sizetf.getValue()).intValue();
				switchSize(size);

			}
		});

		sizesl.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {
				int size = sizesl.getValue();
				switchSize(size);
			}
		});

		familiescb.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				Family family2 = (Family) familiescb.getSelectedItem();
				if (family.getStep() != family2.getStep()) {
					int result = JOptionPane.showConfirmDialog(
							ParticlePickerJFrame.this,
							String.format(
									"Selecting family %s will take you to %s mode."
											+ "\nAre you sure you want to continue?",
									family2.getName(), family2.getStep()
											.toString()), "Message",
							JOptionPane.YES_NO_OPTION);
					if (result == JOptionPane.NO_OPTION)
						return;
					setStep(family2.getStep(), false);
				}
				family = family2;
				color = (family.getColor());
				colorbt.setIcon(new ColorIcon(color));
				sizesl.setValue(family.getSize());
				// ParticlePickerJFrame.this.micrographsmd.fireTableStructureChanged();
				ParticlePickerJFrame.this.micrographsmd.fireTableDataChanged();
				ParticlePickerJFrame.this.micrographstb.getColumnModel()
						.getColumn(1).setHeaderValue(family.getName());
				actionsbt.setText(getFamilyData().getAction());
				actionsbt.setVisible(getFamilyData().isActionAvailable());
			}
		});
		actionsbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				if (actionsbt.getText().equals(State.Autopick.toString())) {
					ppicker.classify(family, micrograph);
					ppicker.loadAutomaticParticles(micrograph);
					getFamilyData().setState(State.Correct);
					canvas.repaint();
				}
				else if(actionsbt.getText().equals(State.Correct.toString())) 
				{
					getFamilyData().setState(State.ReadOnly);
					saveChanges();
					ppicker.persistAutomaticParticles(getFamilyData());
					ppicker.correct(family, micrograph);
					int next = ppicker.getNextFreeMicrograph(family);
					if(next != -1 )
						micrographstb.setRowSelectionInterval(next, next);
					else
						actionsbt.setVisible(false);
					
				}
				actionsbt.setText(getFamilyData().getAction());
			}
		});
	}
	
	MicrographFamilyData getFamilyData()
	{
		return micrograph.getFamilyData(family);
	}

	private void setStep(Step step, boolean update) {
		nextbt.setText("Go To " + ParticlePicker.nextStep(step).toString());
		steplb.setText(step.toString());
		if (update) {
			family.goToNextStep();
			setChanged(true);
			canvas.repaint();
		}
		sizesl.setEnabled(step == Step.Manual);
		sizetf.setEnabled(step == Step.Manual);
		editfamiliesmn.setEnabled(step == Step.Manual);

	}
	
	
	class ColorActionListener implements ActionListener {
		JColorChooser colorChooser;

		@Override
		public void actionPerformed(ActionEvent e) {
			// Set up the dialog that the button brings up.
			colorChooser = new JColorChooser();
			JDialog dialog = JColorChooser.createDialog(colorbt,
					"Pick a Color", true, // modal
					colorChooser, new ActionListener() {

						@Override
						public void actionPerformed(ActionEvent e) {
							family.setColor(ColorActionListener.this.colorChooser
									.getColor());
							updateFamilyColor();
						}
					}, // OK button handler
					null); // no CANCEL button handler
			WindowUtils.centerScreen(position, dialog);
			dialog.setVisible(true);
		}
	}

	private void initSymbolPane() {

		symbolpn = new JPanel();
		symbolpn.setBorder(BorderFactory.createTitledBorder("Symbol"));
		ShapeItemListener shapelistener = new ShapeItemListener();

		circlechb = new JCheckBox(Shape.Circle.toString());
		circlechb.getModel().setSelected(true);
		circlechb.addItemListener(shapelistener);

		rectanglechb = new JCheckBox(Shape.Rectangle.toString());
		rectanglechb.getModel().setSelected(true);
		rectanglechb.addItemListener(shapelistener);

		centerchb = new JCheckBox(Shape.Center.toString());
		centerchb.addItemListener(shapelistener);

		// onlylastchb = new JCheckBox(Shape.OnlyLast.toString());
		// onlylastchb.addItemListener(shapelistener);

		symbolpn.add(circlechb);
		symbolpn.add(rectanglechb);
		symbolpn.add(centerchb);
		// symbolpn.add(onlylastchb);
	}

	class ShapeItemListener implements ItemListener {
		@Override
		public void itemStateChanged(ItemEvent e) {
			canvas.repaint();
		}
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		try {
			activemacro = ((JMenuItem) e.getSource()).getText();
			IJ.run(activemacro);
		} catch (Exception ex) {
			ex.printStackTrace();
			JOptionPane.showMessageDialog(this, ex.getMessage());
		}

	}

	private void initMicrographsPane() {
		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(0, 5, 0, 5);
		constraints.anchor = GridBagConstraints.WEST;
		micrographpn = new JPanel(new GridBagLayout());
		micrographpn.setBorder(BorderFactory.createTitledBorder("Micrograph"));
		JScrollPane sp = new JScrollPane();
		JPanel ctfpn = new JPanel();
		ctfpn.setBorder(BorderFactory.createTitledBorder("CTF"));
		iconlb = new JLabel();
		ctfpn.add(iconlb);
		micrographsmd = new MicrographsTableModel(this);
		micrographstb = new JTable(micrographsmd);
		micrographstb.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		micrographstb.getColumnModel().getColumn(0).setPreferredWidth(180);
		micrographstb.getColumnModel().getColumn(1).setPreferredWidth(50);
		micrographstb.getColumnModel().getColumn(2).setPreferredWidth(70);
		micrographstb
				.setPreferredScrollableViewportSize(new Dimension(300, 150));
		micrographstb
				.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		micrographstb.getSelectionModel().addListSelectionListener(
				new ListSelectionListener() {

					@Override
					public void valueChanged(ListSelectionEvent e) {
						if (e.getValueIsAdjusting())
							return;
						int index = ParticlePickerJFrame.this.micrographstb
								.getSelectedRow();
						if (index == -1)
							return;// Probably from fireTableDataChanged raised
									// by me.
						micrograph = (Micrograph) ppicker.getMicrographs().get(
								index);
						
						initializeCanvas();
						ParticlePickerJFrame.this.iconlb.setIcon(micrograph
								.getCTFIcon());
						actionsbt.setText(getFamilyData().getAction());
						actionsbt.setVisible(getFamilyData().isActionAvailable());
					}
				});
		micrographstb.getSelectionModel().setSelectionInterval(0, 0);
		sp.setViewportView(micrographstb);
		micrographpn.add(sp,
				WindowUtils.updateConstraints(constraints, 0, 0, 1));

		micrographpn.add(ctfpn,
				WindowUtils.updateConstraints(constraints, 1, 0, 1));
	}

	void initializeCanvas() {
		Micrograph micrograph = getMicrograph();
		if (iw == null) {
			canvas = new ParticlePickerCanvas(this);
			iw = new ImageWindow(micrograph.getImage(), canvas);
		} else {
			canvas.updateMicrograph();
		}
		iw.setTitle(micrograph.getName());
		canvas.setName(micrograph.getName());
	}

	private void saveChanges() {

		ppicker.persistFamilies();
		ppicker.persistMicrographsData();
		if(family.getStep() == Step.Supervised)
			ppicker.persistAutomaticParticles(getFamilyData());
		setChanged(false);
	}

	void updateFamilyColor() {
		setChanged(true);
		color = family.getColor();
		colorbt.setIcon(new ColorIcon(color));
		canvas.repaint();
	}

	void updateFamilyComboBox() {
		setChanged(true);
		Family item = (Family) familiescb.getSelectedItem();
		DefaultComboBoxModel model = new DefaultComboBoxModel(ppicker
				.getFamilies().toArray());
		familiescb.setModel(model);
		familiescb.setSelectedItem(item);
		pack();
	}

	void switchSize(int size) {
		sizetf.setText(Integer.toString(size));
		sizesl.setValue(size);
		canvas.repaint();
		family.setSize(size);
		setChanged(true);
	}

	public void addGroup(Family g) {
		if (ppicker.existsFamilyName(g.getName()))
			throw new IllegalArgumentException(
					Constants.getAlreadyExistsGroupNameMsg(g.getName()));
		ppicker.getFamilies().add(g);
		updateFamilyComboBox();
	}

	public void removeFamily(Family family) {
		ppicker.removeFamily(family);
		updateFamilyComboBox();
	}

	void setChanged(boolean changed) {
		this.changed = changed;
		savemi.setEnabled(changed);
	}

	void updateMicrographsModel() {
		micrographsmd.fireTableDataChanged();
	}

}
