package xmipp.particlepicker.training.gui;

import ij.gui.ImageWindow;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.text.NumberFormat;
import java.util.List;
import java.util.logging.Level;
import javax.swing.BorderFactory;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JColorChooser;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFormattedTextField;
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
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import Jama.examples.MagicSquareExample;
import xmipp.particlepicker.Family;
import xmipp.particlepicker.Format;
import xmipp.particlepicker.Micrograph;
import xmipp.particlepicker.ParticlePickerCanvas;
import xmipp.particlepicker.ParticlePickerJFrame;
import xmipp.particlepicker.ParticlesJDialog;
import xmipp.particlepicker.tiltpair.gui.TiltPairParticlesJDialog;
import xmipp.particlepicker.training.gui.TrainingCanvas;
import xmipp.particlepicker.training.model.FamilyState;
import xmipp.particlepicker.training.model.MicrographFamilyData;
import xmipp.particlepicker.training.model.MicrographFamilyState;
import xmipp.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.particlepicker.training.model.TrainingMicrograph;
import xmipp.particlepicker.training.model.TrainingParticle;
import xmipp.particlepicker.training.model.TrainingPicker;
import xmipp.utils.ColorIcon;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippResource;
import xmipp.utils.XmippWindowUtil;

public class TrainingPickerJFrame extends ParticlePickerJFrame {

	private TrainingCanvas canvas;
	private JMenuBar mb;
	private JComboBox familiescb;
	private TrainingPicker ppicker;
	private JPanel familypn;
	private JPanel micrographpn;
	private MicrographsTableModel micrographsmd;
	private TrainingMicrograph micrograph;
	private float positionx;
	private JLabel iconlb;
	private JLabel steplb;
	private JButton actionsbt;
	private JMenuItem editfamiliesmi;
	private int index;
	private JLabel manuallb;
	private JLabel autolb;
	private JSlider thresholdsl;
	private JPanel thresholdpn;
	private JFormattedTextField thresholdtf;
	private Family family;
	private ImageWindow iw;
	private JMenuItem templatesmi;
	TemplatesJDialog templatesdialog;

	@Override
	public TrainingPicker getParticlePicker() {
		return ppicker;
	}

	public TrainingPickerJFrame(TrainingPicker picker) {

		super(picker);
		try {
			this.ppicker = picker;
			initComponents();
			if (ppicker.getMode() == FamilyState.ReadOnly)
				enableEdition(false);
			else if (family.getStep() == FamilyState.Manual
					&& ppicker.getMode() == FamilyState.Supervised)
				goToNextStep();
		} catch (IllegalArgumentException ex) {
			close();
			throw ex;
		}
	}

	public TrainingMicrograph getMicrograph() {
		return micrograph;
	}

	private void initComponents() {
		try {
			setResizable(false);
			setTitle("Xmipp Particle Picker - " + ppicker.getMode());
			initMenuBar();
			setJMenuBar(mb);

			GridBagConstraints constraints = new GridBagConstraints();
			constraints.insets = new Insets(0, 5, 0, 5);
			constraints.anchor = GridBagConstraints.WEST;
			setLayout(new GridBagLayout());

			initImagePane();
			add(imagepn, XmippWindowUtil.getConstraints(constraints, 0, 1, 3));
			
			initFamilyPane();
			add(familypn, XmippWindowUtil.getConstraints(constraints, 0, 2, 3));

			

			initMicrographsPane();
			add(micrographpn,
					XmippWindowUtil.getConstraints(constraints, 0, 3, 3));

			pack();
			positionx = 0.995f;
			XmippWindowUtil.setLocation(positionx, 0.25f, this);
			setVisible(true);
		} catch (Exception e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	protected void enableEdition(boolean enable) {
		super.enableEdition(enable);

		editfamiliesmi.setEnabled(enable);
	}

	public void initMenuBar() {
		mb = new JMenuBar();

		// Setting menus

		JMenu windowmn = new JMenu("Window");
		JMenu helpmn = new JMenu("Help");
		mb.add(filemn);
		mb.add(filtersmn);
		mb.add(windowmn);
		mb.add(helpmn);
		//importffilemi.setText("Import from File...");

		windowmn.add(pmi);
		windowmn.add(ijmi);

		templatesmi = new JMenuItem("Templates");
		editfamiliesmi = new JMenuItem("Edit Families",
				XmippResource.getIcon("edit.gif"));
		windowmn.add(editfamiliesmi);
		windowmn.add(templatesmi);
		helpmn.add(hcontentsmi);

		// Setting menu item listeners

		editfamiliesmi.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				new EditFamiliesJDialog(TrainingPickerJFrame.this, true);

			}
		});
		templatesmi.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				loadTemplates();

			}
		});

	}

	public void loadTemplates() {
		try {
			if (templatesdialog == null)
				templatesdialog = new TemplatesJDialog(
						TrainingPickerJFrame.this);
			else {

				templatesdialog.loadTemplates(true);
				templatesdialog.setVisible(true);
			}
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this, ex.getMessage());
			if (templatesdialog != null)
				templatesdialog.close();
			templatesdialog = null;
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
		family = ppicker.getFamily();
		familiescb.setSelectedItem(family);
		familiescb.setEnabled(ppicker.getMode() == FamilyState.Manual
				|| ppicker.getMode() == FamilyState.ReadOnly);// available
																// families are
																// marked by
																// pickers at
																// start

		if (ppicker.getMode() == FamilyState.Manual
				&& family.getStep() != FamilyState.Manual)
			throw new IllegalArgumentException(
					String.format(
							"Application not enabled for %s mode. Family %s could not be loaded",
							family.getStep(), family.getName()));

		fieldspn.add(familiescb);

		// Setting color
		initColorPane();
		fieldspn.add(colorpn);

		// Setting slider
		initSizePane();
		fieldspn.add(sizepn);

		familypn.add(fieldspn, 0);
		JPanel steppn = new JPanel(new FlowLayout(FlowLayout.LEFT));
		steppn.add(new JLabel("Step:"));
		FamilyState step = family.getStep();
		steplb = new JLabel();
		steppn.add(steplb);

		index = ppicker.getNextFreeMicrograph(0);
		if (index == -1)
			index = 0;
		micrograph = ppicker.getMicrographs().get(index);

		initThresholdPane();
		steppn.add(thresholdpn);
		actionsbt = XmippWindowUtil.getTextButton("", new ActionListener() {

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
		setStep(step);
		steppn.add(actionsbt);

		colorbt.addActionListener(new ColorActionListener());

		familypn.add(steppn, 1);

		familiescb.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				Family family2 = (Family) familiescb.getSelectedItem();
				if (family == family2)
					return;
				// You can only switch between different states for readonly
				// mode. Besides switching will be availabe on manual and
				// readonly modes
				if (family.getStep() != family2.getStep()
						&& ppicker.getMode() != FamilyState.ReadOnly) {
					familiescb.setSelectedItem(family);
					JOptionPane
							.showMessageDialog(
									TrainingPickerJFrame.this,
									String.format(
											"Application not enabled for %s mode. Family %s could not be loaded",
											FamilyState.Supervised,
											family2.getName()));
					return;

				}
				family = family2;
				ppicker.setFamily(family);
				index = ppicker.getNextFreeMicrograph(0);
				micrographstb.getSelectionModel().setSelectionInterval(index,
						index);
				color = (family.getColor());
				colorbt.setIcon(new ColorIcon(color));
				sizesl.setValue(family.getSize());
				updateMicrographsModel();
				micrographstb.getColumnModel()
						.getColumn(micrographsmd.getParticlesPosition())
						.setHeaderValue(family.getName());
				setThresholdValue(0);
			}
		});

	}

	protected void goToNextStep() {

		family.validateNextStep(ppicker);// throws error

		if (family.getStep() == FamilyState.Manual)
			train();

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
			XmippWindowUtil.setLocation(positionx, 0.25f, dialog);
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
				setThresholdValue(threshold);

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

	private void setThresholdValue(double threshold) {
		if (Math.abs(threshold) <= 100) {
			thresholdsl.setValue((int) threshold);
			setThresholdChanges();
		}
	}

	private void setThresholdChanges() {
		setChanged(true);
		updateMicrographsModel();
		canvas.repaint();
		actionsbt.setVisible(getFamilyData().isActionVisible(getThreshold()));
		if (particlesdialog != null)
			loadParticles();

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
				TitledBorder.CENTER, TitledBorder.BELOW_BOTTOM));
		iconlb = new JLabel();
		ctfpn.add(iconlb);
		micrographsmd = new MicrographsTableModel(this);
		micrographstb.setModel(micrographsmd);
		formatMicrographsTable();

		sp.setViewportView(micrographstb);
		micrographpn.add(sp,
				XmippWindowUtil.getConstraints(constraints, 0, 0, 1));
		micrographpn.add(ctfpn,
				XmippWindowUtil.getConstraints(constraints, 1, 0, 1));
		JPanel infopn = new JPanel();
		manuallb = new JLabel(Integer.toString(ppicker
				.getManualParticlesNumber(family)));
		autolb = new JLabel(Integer.toString(ppicker.getAutomaticNumber(family,
				getThreshold())));
		infopn.add(new JLabel("Manual:"));
		infopn.add(manuallb);
		infopn.add(new JLabel("Automatic:"));
		infopn.add(autolb);
		micrographpn.add(infopn,
				XmippWindowUtil.getConstraints(constraints, 0, 1, 1));
		JPanel buttonspn = new JPanel(new FlowLayout(FlowLayout.LEFT));

		buttonspn.add(resetbt);
		micrographpn.add(buttonspn,
				XmippWindowUtil.getConstraints(constraints, 0, 2, 2));

		micrographstb.getSelectionModel().setSelectionInterval(index, index);

	}

	protected void loadMicrograph() {
		if (TrainingPickerJFrame.this.micrographstb.getSelectedRow() == -1)
			return;// Probably from fireTableDataChanged raised
		if (index == TrainingPickerJFrame.this.micrographstb.getSelectedRow()
				&& iw != null && iw.isVisible())// same micrograph open
			return;
		index = TrainingPickerJFrame.this.micrographstb.getSelectedRow();
		// by me.
		micrograph.releaseImage();
		micrograph = (TrainingMicrograph) ppicker.getMicrographs().get(index);

		initializeCanvas();
		TrainingPickerJFrame.this.iconlb.setIcon(micrograph.getCTFIcon());
		actionsbt.setText(getFamilyData().getAction());
		actionsbt.setVisible(getFamilyData().isActionVisible(getThreshold()));
		thresholdpn
				.setVisible(getFamilyData().getState() == MicrographFamilyState.Correct);
		pack();
		ppicker.saveData(getMicrograph());// Saving changes when switching micrographs, by Coss suggestion
		if (particlesdialog != null)
			loadParticles();

	}

	protected void resetMicrograph() {
		ppicker.resetFamilyData(getFamilyData());
		canvas.setActive(null);
		updateMicrographsModel();
		setState(MicrographFamilyState.Available);

	}

	private void setState(MicrographFamilyState state) {
		setChanged(true);
		getFamilyData().setState(state);
		actionsbt.setText(getFamilyData().getAction());
		// if (getFamilyData().getState() == MicrographFamilyState.Correct)
		// actionsbt.setEnabled(false);// enabled only after doing corrections
		saveChanges();// to keep consistence between files of automatic picker
						// and mines
		thresholdpn.setVisible(state == MicrographFamilyState.Correct);
		updateMicrographsModel();
		pack();
	}

	public MicrographFamilyData getFamilyData() {
		MicrographFamilyData mfd = micrograph.getFamilyData(family);
		return mfd;
	}

	private void setStep(FamilyState step) {
		if (micrographsmd != null) {
			updateMicrographsModel();
			canvas.repaint();// paints only current class in supervised mode
		}

		steplb.setText(step.toString());
		sizesl.setEnabled(step == FamilyState.Manual);
		sizetf.setEnabled(step == FamilyState.Manual);
		editfamiliesmi.setEnabled(step == FamilyState.Manual);
		actionsbt.setText(getFamilyData().getAction());
		actionsbt.setVisible(getFamilyData().isActionVisible(getThreshold()));
		thresholdpn
				.setVisible(getFamilyData().getState() == MicrographFamilyState.Correct);
		pack();

	}

	protected void initializeCanvas() {
		
		if (canvas == null) {
			canvas = new TrainingCanvas(this);
			iw = new ImageWindow(micrograph.getImagePlus(ppicker.getFilters()),
					canvas);
			iw.setTitle(micrograph.getName());
		} else {
			canvas.updateMicrograph();
			// seems to keep previous window instead of creating a new one
			iw = new ImageWindow(canvas.getImage(), canvas);

		}
		micrograph.runImageJFilters(ppicker.getFilters());
		
		double zoom = Double.parseDouble(usezoombt.getText());
		if(zoom == -1. || ( zoom != -1. && !usezoombt.isSelected()))//setting canvas magnification
		{
			zoom = canvas.getMagnification();
			usezoombt.setText(String.format("%.2f", zoom));
		}
		else if(usezoombt.isSelected())
			canvas.setZoom(zoom);
	}

	private void formatMicrographsTable() {
		micrographstb.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		micrographstb.getColumnModel().getColumn(0).setPreferredWidth(35);
		micrographstb.getColumnModel().getColumn(1).setPreferredWidth(245);
		micrographstb.getColumnModel().getColumn(2).setPreferredWidth(70);
		micrographstb.getColumnModel().getColumn(3).setPreferredWidth(70);
		micrographstb
				.setPreferredScrollableViewportSize(new Dimension(420, 304));
		micrographstb.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		if (index != -1)
			micrographstb.setRowSelectionInterval(index, index);
	}

	protected void saveChanges() {

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
		micrographsmd.fireTableStructureChanged();

		formatMicrographsTable();
		pack();
	}

	public void addFamily(Family g) {
		if (ppicker.existsFamilyName(g.getName()))
			throw new IllegalArgumentException(
					XmippMessage.getAlreadyExistsGroupNameMsg(g.getName()));
		ppicker.getFamilies().add(g);
		updateFamilyComboBox();
	}

	public void removeFamily(Family family) {
		ppicker.removeFamily(family);
		updateFamilyComboBox();
	}

	public void setChanged(boolean changed) {
		ppicker.setChanged(changed);
		savemi.setEnabled(changed);
	}

	public void updateMicrographsModel() {
		super.updateMicrographsModel();
		if (templatesdialog != null)
			loadTemplates();
		micrographsmd.fireTableRowsUpdated(index, index);
		micrographstb.setRowSelectionInterval(index, index);
		manuallb.setText(Integer.toString(ppicker
				.getManualParticlesNumber(family)));
		autolb.setText(Integer.toString(ppicker.getAutomaticNumber(family,
				getThreshold())));
		actionsbt.setVisible(getFamilyData().isActionVisible(getThreshold()));
	}

	public ParticlePickerCanvas getCanvas() {
		return canvas;
	}

	private void train() {

		family.goToNextStep(ppicker);// validate and change state if posible
		// setChanged(true);
		setStep(FamilyState.Supervised);// change visual appearance
		saveChanges();// persist changes

		try {
			canvas.setEnabled(false);
			XmippWindowUtil.blockGUI(getRootPane(), "Training...");
			Thread t = new Thread(new Runnable() {

				public void run() {
					MicrographFamilyData mfd;
					String args;
					SupervisedParticlePicker sppicker = ((SupervisedParticlePicker) ppicker);
					for (TrainingMicrograph micrograph : ppicker
							.getMicrographs()) {
						mfd = micrograph.getFamilyData(family);
						if (!mfd.isEmpty()) {
							args = sppicker
									.getBuildInvariantCommandLineArgs(mfd);
							ppicker.runXmippProgram(
									"xmipp_micrograph_automatic_picking", args);
							System.out.println(args);
						}
					}

					args = sppicker.getTrainCommandLineArgs();
					System.out.println(args);
					ppicker.runXmippProgram(
							"xmipp_micrograph_automatic_picking", args);
					int next = ppicker.getNextFreeMicrograph(index);
					if (next != -1)
						micrographstb.setRowSelectionInterval(next, next);
					XmippWindowUtil.releaseGUI(getRootPane());
					canvas.setEnabled(true);
				}
			});
			t.start();

		} catch (Exception e) {
			TrainingPicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	private void autopick() {
		setState(MicrographFamilyState.Autopick);

		final String fargs = ((SupervisedParticlePicker) ppicker)
				.getAutopickCommandLineArgs(getFamilyData());
		System.out.println(fargs);
		try {
			canvas.setEnabled(false);
			XmippWindowUtil.blockGUI(getRootPane(), "Autopicking...");
			Thread t = new Thread(new Runnable() {

				public void run() {
					ppicker.runXmippProgram(
							"xmipp_micrograph_automatic_picking", fargs);
					ppicker.loadAutomaticParticles(getFamilyData());
					setState(MicrographFamilyState.Correct);
					canvas.repaint();
					canvas.setEnabled(true);
					XmippWindowUtil.releaseGUI(getRootPane());
				}
			});
			t.start();

		} catch (Exception e) {
			TrainingPicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	private void correct() {
		getFamilyData().deleteBelowThreshold(getThreshold());
		setState(MicrographFamilyState.ReadOnly);
		ppicker.persistAutomaticParticles(getFamilyData());

		try {
			canvas.setEnabled(false);
			XmippWindowUtil.blockGUI(getRootPane(), "Correcting...");

			Thread t = new Thread(new Runnable() {
				public void run() {
					SupervisedParticlePicker sppicker = ((SupervisedParticlePicker) ppicker);
					String args;
					args = sppicker
							.getBuildInvariantCommandLineArgs(getFamilyData());
					sppicker.runXmippProgram(
							"xmipp_micrograph_automatic_picking", args);// build
																		// invariants
					args = sppicker.getCorrectCommandLineArgs(getFamilyData());
					ppicker.runXmippProgram(
							"xmipp_micrograph_automatic_picking", args);// correct
					int next = ppicker.getNextFreeMicrograph(index + 1);
					if (next != -1)
						micrographstb.setRowSelectionInterval(next, next);
					else
						actionsbt.setVisible(false);
					canvas.setEnabled(true);
					XmippWindowUtil.releaseGUI(getRootPane());
				}
			});
			t.start();

		} catch (Exception e) {
			TrainingPicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public double getThreshold() {
		return thresholdsl.getValue() / 100.0;
	}

	@Override
	public List<? extends TrainingParticle> getAvailableParticles() {
		return getFamilyData().getAvailableParticles(getThreshold());
	}

	@Override
	public boolean isPickingAvailable(MouseEvent e) {
		if (!super.isPickingAvailable(e))
			return false;
		return getFamilyData().isPickingAvailable();
	}

	@Override
	public void changeShapes() {
		canvas.repaint();

	}

	@Override
	protected void reloadImage() {
		getCanvas().getMicrograph().releaseImage();
		getCanvas().updateMicrographData();

	}

	
	public void importParticlesFromFile(Format format, String file) {
		String filename = Micrograph.getName(file, 1);
		if(!filename.equals(getMicrograph().getName()))//validating you want use this file for this micrograph with different name
		{
			String msg = String.format("Are you sure you want to import data from file\n%s to micrograph %s ?", file, getMicrograph().getName());
			int result = JOptionPane.showConfirmDialog(this, msg);
			if(result != JOptionPane.YES_OPTION)
				return;
		}
		MicrographFamilyData mfd = getFamilyData();
		mfd.reset();
		ppicker.importParticlesFromFile(file, format, mfd.getMicrograph());
		setChanged(true);
		getCanvas().repaint();
		updateMicrographsModel();
		updateSize(family.getSize());
		canvas.setActive(null);
	}

	@Override
	public boolean isValidSize(int size) {
		for (TrainingParticle p : getFamilyData().getParticles())
			if (!micrograph.fits(p.getX(), p.getY(), size))
				return false;
		return true;
	}

	@Override
	protected void openHelpURl() {
		XmippWindowUtil
				.openURI("http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Micrograph_particle_picking_v3");

	}

	public void updateTemplates() {
		ppicker.updateFamilyTemplates(family);

	}

	public void updateSize(int size) {
		super.updateSize(size);
		if (templatesdialog != null)
			loadTemplates();
	}
	
	@Override
	protected void resetData(){
		getFamilyData().reset();
	}
}
