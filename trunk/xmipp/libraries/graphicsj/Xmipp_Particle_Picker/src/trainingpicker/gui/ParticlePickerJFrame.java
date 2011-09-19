package trainingpicker.gui;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.ImageWindow;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;

import javax.swing.BorderFactory;
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
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import trainingpicker.model.Constants;
import trainingpicker.model.Family;
import trainingpicker.model.FamilyState;
import trainingpicker.model.TrainingMicrograph;
import trainingpicker.model.MicrographFamilyData;
import trainingpicker.model.MicrographFamilyState;
import trainingpicker.model.TrainingParticle;
import trainingpicker.model.ParticlePicker;
import trainingpicker.model.SupervisedParticlePicker;
import trainingpicker.model.XmippJ;

import xmipp.Program;


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

	private JMenuItem savemi;
	private MicrographsTableModel micrographsmd;
	TrainingMicrograph micrograph;
	private JButton nextbt;
	private JButton colorbt;
	private double position;
	private JLabel iconlb;
	private JLabel steplb;
	private JButton actionsbt;
	private Family family;
	private JMenuItem editfamiliesmn;
	private int index;
	private JButton resetbt;
	private JLabel manuallb;
	private JLabel autolb;
	private JSlider thresholdsl;
	private JPanel thresholdpn;
	private JFormattedTextField thresholdtf;
	private String tool = "Particle Picker Tool";

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

		if (IJ.getInstance() == null)
			return Tool.PICKER;
		if (IJ.getToolName().equalsIgnoreCase(tool))
			return Tool.PICKER;
		return Tool.IMAGEJ;
	}

	public ParticlePicker getParticlePicker() {
		return ppicker;
	}

	public Family getFamily() {
		return family;
	}

	public ParticlePickerJFrame(ParticlePicker ppicker) {

		this.ppicker = ppicker;
		initComponents();
	}

	public TrainingMicrograph getMicrograph() {
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
				if (ppicker.isChanged()) {
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
		setTitle("Xmipp Particle Picker - " + ppicker.getMode());
		initMenuBar();
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

	public void initMenuBar() {
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
		savemi.setEnabled(ppicker.isChanged());
		filemn.add(savemi);
		savemi.setMnemonic('S');
		savemi.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S,
				InputEvent.CTRL_DOWN_MASK));

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

					new ImageJ();
					IJ.run("Install...",
							"install="
									+ ParticlePicker
											.getXmippPath("external/imagej/macros/ParticlePicker.txt"));
					IJ.setTool(tool);
				}
				// IJ.getInstance().setVisible(true);
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
				for (TrainingParticle p : getMicrograph().getFamilyData(family)
						.getManualParticles())
					imgs.add(p.getImage());
				String filename = XmippJ.saveTempImageStack(imgs);
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
		hcontentsmi.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				try {
					WindowUtils
							.openURI("http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParticlePicker");
				} catch (Exception ex) {
					JOptionPane.showMessageDialog(ParticlePickerJFrame.this,
							ex.getMessage());
				}
			}
		});

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

	private void initFamilyPane() {
		familypn = new JPanel();
		GridLayout gl = new GridLayout(2, 1);
		familypn.setLayout(gl);

		familypn.setBorder(BorderFactory.createTitledBorder("Family"));

		JPanel fieldspn = new JPanel(new FlowLayout(FlowLayout.LEFT));

		// Setting combo
		fieldspn.add(new JLabel("Name:"));
		familiescb = new JComboBox(ppicker.getFamilies().toArray());
		familiescb.setEnabled(ppicker.getMode() != FamilyState.Review);

		family = (Family) familiescb.getSelectedItem();
		if (ppicker.getMode() == FamilyState.Manual
				&& family.getStep() == FamilyState.Supervised)
			throw new IllegalArgumentException(
					String.format(
							"Application not enabled for %s mode. Family %s could not be loaded",
							FamilyState.Supervised, family.getName()));

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
		sizesl.setMajorTickSpacing(250);
		sizesl.setPaintTicks(true);
		sizesl.setPaintLabels(true);
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
		FamilyState step = family.getStep();
		steplb = new JLabel();
		steppn.add(steplb);

		nextbt = new JButton();
		nextbt.setVisible(ppicker.getMode() == FamilyState.Supervised);
		steppn.add(nextbt);

		index = ppicker.getNextFreeMicrograph(family);
		if (index == -1)
			index = 0;
		micrograph = ppicker.getMicrographs().get(index);

		initThresholdPane();
		steppn.add(thresholdpn);
		actionsbt = new JButton();
		setStep(step);
		steppn.add(actionsbt);

		colorbt.addActionListener(new ColorActionListener());
		nextbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				try {
					family.validateNextStep(ppicker);
				} catch (Exception ex) {
					JOptionPane.showMessageDialog(ParticlePickerJFrame.this,
							ex.getMessage());
					return;
				}
				if (family.getStep() == FamilyState.Manual) {
					int result = JOptionPane
							.showConfirmDialog(
									ParticlePickerJFrame.this,
									"Data provided from manual picking will be used to train software and start Supervised mode."
											+ "\nManual mode will be disabled. Are you sure you want to continue?",
									"Message", JOptionPane.YES_NO_OPTION);

					if (result == JOptionPane.NO_OPTION)
						return;
					train();

				} else if (family.getStep() == FamilyState.Supervised) {
					int result = JOptionPane.showConfirmDialog(
							ParticlePickerJFrame.this,
							"Model generated during Supervised mode will be dismissed."
									+ "\nAre you sure you want to continue?",
							"Message", JOptionPane.YES_NO_OPTION);

					if (result == JOptionPane.NO_OPTION)
						return;
					family.goToPreviousStep();
					((SupervisedParticlePicker)ppicker).resetModel(family);
					setStep(FamilyState.Manual);

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
				if (ppicker.getMode() == FamilyState.Manual
						&& family2.getStep() == FamilyState.Supervised) {
					familiescb.setSelectedItem(family);
					JOptionPane
							.showMessageDialog(
									ParticlePickerJFrame.this,
									String.format(
											"Application not enabled for %s mode. Family %s could not be loaded",
											FamilyState.Supervised,
											family2.getName()));
					return;
				}
				if (family.getStep() != family2.getStep()) {
					int result = JOptionPane.showConfirmDialog(
							ParticlePickerJFrame.this,
							String.format(
									"Selecting family %s will take you to %s mode."
											+ "\nAre you sure you want to continue?",
									family2.getName(), family2.getStep()
											.toString()), "Message",
							JOptionPane.YES_NO_OPTION);
					if (result == JOptionPane.NO_OPTION) {
						familiescb.setSelectedItem(family);
						return;
					}
					setStep(family2.getStep());
				}
				family = family2;
				color = (family.getColor());
				colorbt.setIcon(new ColorIcon(color));
				sizesl.setValue(family.getSize());
				updateMicrographsModel();
				micrographstb.getColumnModel()
						.getColumn(micrographsmd.getParticlesPosition())
						.setHeaderValue(family.getName());
			}
		});
		actionsbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				if (actionsbt.getText().equals(
						MicrographFamilyState.Autopick.toString()))
					autopick();
				else if (actionsbt.getText().equals(
						MicrographFamilyState.Correct.toString()))
					correct();
			}
		});

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

	private void initThresholdPane() {
		thresholdpn = new JPanel();
		thresholdpn.add(new JLabel("Threshold:"));
		thresholdsl = new JSlider(0, 100, 0);
		thresholdsl.setPaintTicks(true);
		thresholdsl.setMajorTickSpacing(50);
		java.util.Hashtable<Integer, JLabel> labelTable = new java.util.Hashtable<Integer, JLabel>();
		labelTable.put(new Integer(100), new JLabel("1"));
		labelTable.put(new Integer(50), new JLabel("0.5"));
		labelTable.put(new Integer(0), new JLabel("0"));
		thresholdsl.setLabelTable(labelTable);
		thresholdsl.setPaintLabels(true);
		int height = (int) thresholdsl.getPreferredSize().getHeight();

		thresholdsl.setPreferredSize(new Dimension(100, height));
		thresholdpn.add(thresholdsl);

		thresholdtf = new JFormattedTextField(NumberFormat.getNumberInstance());
		;
		thresholdtf.setColumns(3);
		thresholdtf.setText(Integer.toString(0));
		thresholdpn.add(thresholdtf);
		thresholdpn
				.setVisible(getFamilyData().getState() == MicrographFamilyState.Correct);
		thresholdtf.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				double threshold = Double.parseDouble(thresholdtf.getText()) * 100;
				if (Math.abs(threshold) <= 100) {
					thresholdsl.setValue((int) threshold);
					setThresholdChanges();
				}

			}
		});

		thresholdsl.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {
				double threshold = (double) thresholdsl.getValue() / 100;
				thresholdtf.setText(String.format("%.2f", threshold));
				setThresholdChanges();
			}
		});
	}

	private void setThresholdChanges() {
		updateMicrographsModel();
		canvas.repaint();
		actionsbt.setVisible(getFamilyData().isActionAvailable(getThreshold()));
	}

	private void initSymbolPane() {

		symbolpn = new JPanel();
		symbolpn.setBorder(BorderFactory.createTitledBorder("Symbol"));
		ShapeItemListener shapelistener = new ShapeItemListener();

		circlechb = new JCheckBox(Shape.Circle.toString());
		circlechb.setSelected(true);
		circlechb.addItemListener(shapelistener);

		rectanglechb = new JCheckBox(Shape.Rectangle.toString());
		rectanglechb.setSelected(true);
		rectanglechb.addItemListener(shapelistener);

		centerchb = new JCheckBox(Shape.Center.toString());
		centerchb.setSelected(true);
		centerchb.addItemListener(shapelistener);

		symbolpn.add(circlechb);
		symbolpn.add(rectanglechb);
		symbolpn.add(centerchb);
	}

	class ShapeItemListener implements ItemListener {
		@Override
		public void itemStateChanged(ItemEvent e) {
			canvas.repaint();
		}
	}

	private void initMicrographsPane() {
		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(0, 5, 0, 5);
		constraints.anchor = GridBagConstraints.NORTHWEST;
		micrographpn = new JPanel(new GridBagLayout());
		micrographpn.setBorder(BorderFactory.createTitledBorder("Micrograph"));
		JScrollPane sp = new JScrollPane();
		JPanel ctfpn = new JPanel();
		ctfpn.setBorder(BorderFactory.createTitledBorder(null, "CTF",
				javax.swing.border.TitledBorder.CENTER,
				javax.swing.border.TitledBorder.BELOW_BOTTOM));
		iconlb = new JLabel();
		ctfpn.add(iconlb);
		micrographsmd = new MicrographsTableModel(this);
		micrographstb = new JTable(micrographsmd);
		micrographstb.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		micrographstb.getColumnModel().getColumn(0).setPreferredWidth(35);
		micrographstb.getColumnModel().getColumn(1).setPreferredWidth(245);
		micrographstb.getColumnModel().getColumn(2).setPreferredWidth(70);
		micrographstb.getColumnModel().getColumn(3).setPreferredWidth(70);
		micrographstb
				.setPreferredScrollableViewportSize(new Dimension(420, 304));
		micrographstb.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

		sp.setViewportView(micrographstb);
		micrographpn.add(sp,
				WindowUtils.updateConstraints(constraints, 0, 0, 1));
		micrographpn.add(ctfpn,
				WindowUtils.updateConstraints(constraints, 1, 0, 1));
		JPanel infopn = new JPanel();
		manuallb = new JLabel(Integer.toString(family.getManualNumber()));
		autolb = new JLabel(Integer.toString(ppicker.getAutomaticNumber(family,
				getThreshold())));
		infopn.add(new JLabel("Manual:"));
		infopn.add(manuallb);
		infopn.add(new JLabel("Automatic:"));
		infopn.add(autolb);
		micrographpn.add(infopn,
				WindowUtils.updateConstraints(constraints, 0, 1, 1));
		JPanel buttonspn = new JPanel(new FlowLayout(FlowLayout.LEFT));
		resetbt = new JButton("Reset");
		buttonspn.add(resetbt);
		micrographpn.add(buttonspn,
				WindowUtils.updateConstraints(constraints, 0, 2, 2));
		resetbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				ppicker.resetFamilyData(getFamilyData());
				setState(MicrographFamilyState.Available);
				canvas.repaint();
				updateMicrographsModel();
				setChanged(true);
			}
		});
		micrographstb.getSelectionModel().addListSelectionListener(
				new ListSelectionListener() {

					@Override
					public void valueChanged(ListSelectionEvent e) {
						if (e.getValueIsAdjusting())
							return;
						if (ParticlePickerJFrame.this.micrographstb
								.getSelectedRow() == -1)
							return;// Probably from fireTableDataChanged raised
						index = ParticlePickerJFrame.this.micrographstb
								.getSelectedRow();
						// by me.
						micrograph.releaseImage();
						micrograph = (TrainingMicrograph) ppicker.getMicrographs().get(
								index);

						initializeCanvas();
						ParticlePickerJFrame.this.iconlb.setIcon(micrograph
								.getCTFIcon());
						actionsbt.setText(getFamilyData().getAction());
						actionsbt.setVisible(getFamilyData().isActionAvailable(
								getThreshold()));
						thresholdpn
								.setVisible(getFamilyData().getState() == MicrographFamilyState.Correct);
						pack();
						saveChanges();//Saving changes when switching micrographs, by Coss suggestion
					}
				});
		micrographstb.getSelectionModel().setSelectionInterval(index, index);

	}

	private void setState(MicrographFamilyState state) {
		getFamilyData().setState(state);
		actionsbt.setText(getFamilyData().getAction());
		saveChanges();// to keep consistence between files of automatic picker
						// and mines
		thresholdpn.setVisible(state == MicrographFamilyState.Correct);
		updateMicrographsModel();
		pack();
	}

	MicrographFamilyData getFamilyData() {
		return micrograph.getFamilyData(family);
	}

	private void setStep(FamilyState step) {
		if (micrographsmd != null) {
			updateMicrographsModel();
			canvas.repaint();// paints only current class in supervised mode
		}
		if (step == FamilyState.Manual)
			nextbt.setText("Go To " + ParticlePicker.nextStep(step).toString());
		else if (step == FamilyState.Supervised)
			nextbt.setText("Go To "
					+ ParticlePicker.previousStep(step).toString());
		else
			nextbt.setVisible(false);
		steplb.setText(step.toString());
		sizesl.setEnabled(step == FamilyState.Manual);
		sizetf.setEnabled(step == FamilyState.Manual);
		editfamiliesmn.setEnabled(step == FamilyState.Manual);
		actionsbt.setText(getFamilyData().getAction());
		actionsbt.setVisible(getFamilyData().isActionAvailable(getThreshold()));
		thresholdpn
				.setVisible(getFamilyData().getState() == MicrographFamilyState.Correct);
		pack();
	}

	void initializeCanvas() {
		if (iw == null) {
			canvas = new ParticlePickerCanvas(this);
			iw = new ImageWindow(micrograph.getImage(), canvas);
		} else {
			canvas.updateMicrograph();
			micrograph.getImage().setWindow(iw);
		}
		iw.setTitle(micrograph.getName());
		canvas.setName(micrograph.getName());
	}

	private void saveChanges() {
		ppicker.saveData();
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

	public void addFamily(Family g) {
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
		ppicker.setChanged(changed);
		savemi.setEnabled(changed);
	}

	void updateMicrographsModel() {
		micrographsmd.fireTableRowsUpdated(index, index);
		micrographstb.setRowSelectionInterval(index, index);
		manuallb.setText(Integer.toString(family.getManualNumber()));
		autolb.setText(Integer.toString(ppicker.getAutomaticNumber(family,
				getThreshold())));
		actionsbt.setVisible(getFamilyData().isActionAvailable(getThreshold()));
	}

	public ImageWindow getImageWindow() {
		return iw;
	}

	public ParticlePickerCanvas getCanvas() {
		return canvas;
	}

	private void train() {
		
		saveChanges();
		family.goToNextStep(ppicker);
		setChanged(true);
		setStep(FamilyState.Supervised);

		String args;
		for (TrainingMicrograph micrograph : ppicker.getMicrographs()) {
			if (!micrograph.getFamilyData(family).isEmpty()) {

				args = String
						.format("-i %s --particleSize %s --model %s --outputRoot %s --mode train %s",
								micrograph.getFile(),// -i
								family.getSize(), // --particleSize
								ppicker.getOutputPath(family.getName()),// --model
								ppicker.getOutputPath(micrograph.getName()), // --outputRoot
								family.getName() + "@"
										+ micrograph.getOFilename());// train
				// parameter
				if (((SupervisedParticlePicker)ppicker).isFastMode())
					args += " --fast";
				if (((SupervisedParticlePicker)ppicker).isIncore())
					args += " --in_core";
				final String fargs = args;
				try {
					final InfiniteProgressPanel glassPane = new InfiniteProgressPanel(
							"training picker...");
					final Component previousGlassPane = ParticlePickerJFrame.this
							.getRootPane().getGlassPane();
					canvas.setEnabled(false);
					ParticlePickerJFrame.this.getRootPane().setGlassPane(
							glassPane);
					glassPane.start();

					Thread t = new Thread(new Runnable() {

						public void run() {

							try {

								Program.runByName(
										"xmipp_micrograph_automatic_picking",
										fargs);
							} catch (Exception e) {
								ParticlePicker.getLogger().log(Level.SEVERE,
										e.getMessage(), e);
								throw new IllegalArgumentException(e);
							}
							glassPane.stop();
							ParticlePickerJFrame.this.getRootPane()
									.setGlassPane(previousGlassPane);
							canvas.setEnabled(true);
							int next = ppicker.getNextFreeMicrograph(family);
							if (next != -1)
								micrographstb.setRowSelectionInterval(next,
										next);

						}
					});
					t.start();

				} catch (Exception e) {
					ParticlePicker.getLogger().log(Level.SEVERE,
							e.getMessage(), e);
					throw new IllegalArgumentException(e.getMessage());
				}
			}
		}
	}

	private void autopick() {
		setState(MicrographFamilyState.Autopick);
		String args;
		args = String
				.format("-i %s --particleSize %s --model %s --outputRoot %s --mode try --thr %s",
						micrograph.getFile(),// -i
						family.getSize(), // --particleSize
						ppicker.getOutputPath(family.getName()),// --model
						ppicker.getOutputPath(micrograph.getName()),// --outputRoot
						((SupervisedParticlePicker)ppicker).getThreads()// --thr
				);

		if (((SupervisedParticlePicker)ppicker).isFastMode())
			args += " --fast";
		if (((SupervisedParticlePicker)ppicker).isIncore())
			args += " --in_core";
		final String fargs = args;
		try {
			final InfiniteProgressPanel glassPane = new InfiniteProgressPanel(
					"autopicking...");
			final Component previousGlassPane = ParticlePickerJFrame.this
					.getRootPane().getGlassPane();
			canvas.setEnabled(false);
			ParticlePickerJFrame.this.getRootPane().setGlassPane(glassPane);
			glassPane.start();

			Thread t = new Thread(new Runnable() {

				public void run() {

					try {

						Program.runByName("xmipp_micrograph_automatic_picking",
								fargs);
					} catch (Exception e) {
						ParticlePicker.getLogger().log(Level.SEVERE,
								e.getMessage(), e);
						throw new IllegalArgumentException(e);
					}
					glassPane.stop();
					ParticlePickerJFrame.this.getRootPane().setGlassPane(
							previousGlassPane);
					ppicker.loadAutomaticParticles(micrograph);
					setState(MicrographFamilyState.Correct);
					canvas.repaint();
					canvas.setEnabled(true);

				}
			});
			t.start();

		} catch (Exception e) {
			ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	private void correct() {
		getFamilyData().deleteBelowThreshold(getThreshold());
		setState(MicrographFamilyState.ReadOnly);
		ppicker.persistAutomaticParticles(getFamilyData());

		String args = String
				.format("-i %s --particleSize %s --model %s --outputRoot %s --mode train ",
						micrograph.getFile(),// -i
						family.getSize(), // --particleSize
						ppicker.getOutputPath(family.getName()),// --model
						ppicker.getOutputPath(micrograph.getName())// --outputRoot
				);

		if (micrograph.getFamilyData(family).getManualParticles().size() > 0)
			args += family.getName() + "@" + micrograph.getOFilename();
		if (((SupervisedParticlePicker)ppicker).isFastMode())
			args += " --fast";
		if (((SupervisedParticlePicker)ppicker).isIncore())
			args += " --in_core";
		final String fargs = args;
		try {
			final InfiniteProgressPanel glassPane = new InfiniteProgressPanel(
					"correcting...");
			final Component previousGlassPane = ParticlePickerJFrame.this
					.getRootPane().getGlassPane();
			canvas.setEnabled(false);
			ParticlePickerJFrame.this.getRootPane().setGlassPane(glassPane);
			glassPane.start();
			Thread t = new Thread(new Runnable() {
				public void run() {
					try {
						Program.runByName("xmipp_micrograph_automatic_picking",
								fargs);
					} catch (Exception e) {
						ParticlePicker.getLogger().log(Level.SEVERE,
								e.getMessage(), e);
						throw new IllegalArgumentException(e);
					}

					glassPane.stop();

					ParticlePickerJFrame.this.getRootPane().setGlassPane(
							previousGlassPane);
					int next = ppicker.getNextFreeMicrograph(family);
					if (next != -1)
						micrographstb.setRowSelectionInterval(next, next);
					else
						actionsbt.setVisible(false);
					canvas.setEnabled(true);

				}
			});
			t.start();

		} catch (Exception e) {
			ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public double getThreshold() {
		return thresholdsl.getValue() / 100.0;
	}

}
