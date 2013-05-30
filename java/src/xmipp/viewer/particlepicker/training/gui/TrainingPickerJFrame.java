package xmipp.viewer.particlepicker.training.gui;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.text.NumberFormat;
import java.util.List;
import java.util.logging.Level;
import javax.swing.BorderFactory;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
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
import xmipp.utils.ColorIcon;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippResource;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.ctf.CTFAnalyzerJFrame;
import xmipp.viewer.particlepicker.Family;
import xmipp.viewer.particlepicker.Format;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.ParticleToTemplatesTask;
import xmipp.viewer.particlepicker.ParticlesJDialog;
import xmipp.viewer.particlepicker.UpdateTemplatesTask;
import xmipp.viewer.particlepicker.training.model.FamilyState;
import xmipp.viewer.particlepicker.training.model.ManualParticlePicker;
import xmipp.viewer.particlepicker.training.model.MicrographFamilyData;
import xmipp.viewer.particlepicker.training.model.MicrographFamilyState;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.viewer.particlepicker.training.model.TrainingMicrograph;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;
import xmipp.viewer.particlepicker.training.model.TrainingPicker;

public class TrainingPickerJFrame extends ParticlePickerJFrame
{

	private TrainingCanvas canvas;
	private JMenuBar mb;
	private JComboBox familiescb;
	private TrainingPicker ppicker;
	private JPanel familypn;
	private JPanel micrographpn;
	private MicrographsTableModel micrographsmd;

	private float positionx;
	private JButton iconbt;
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
	private JFormattedTextField autopickpercenttf;
	private JPanel autopickpercentpn;

	private JMenuItem templatesmi;
	TemplatesJDialog templatesdialog;
	private JCheckBox centerpickchb;
	private JButton editfamiliesbt;

	@Override
	public TrainingPicker getParticlePicker()
	{
		return ppicker;
	}

	public TrainingPickerJFrame(TrainingPicker picker)
	{

		super(picker);
		try
		{
			this.ppicker = picker;
			initComponents();
			if (ppicker.getMode() == FamilyState.ReadOnly)
				enableEdition(false);
			else if (family.getStep() == FamilyState.Manual && ppicker.getMode() == FamilyState.Supervised)
				goToNextStep();
		}
		catch (IllegalArgumentException ex)
		{
			close();
			throw ex;
		}
	}

	public TrainingMicrograph getMicrograph()
	{
		return ppicker.getMicrograph();
	}

	private void initComponents()
	{
		try
		{
			setResizable(false);
			setTitle("Xmipp Particle Picker - " + ppicker.getMode());
			initMenuBar();
			setJMenuBar(mb);

			GridBagConstraints constraints = new GridBagConstraints();
			constraints.insets = new Insets(0, 5, 0, 5);
			constraints.anchor = GridBagConstraints.WEST;
			setLayout(new GridBagLayout());

			initImagePane();
			add(imagepn, XmippWindowUtil.getConstraints(constraints, 0, 1));

			initFamilyPane();
			add(familypn, XmippWindowUtil.getConstraints(constraints, 0, 2));

			initMicrographsPane();
			add(micrographpn, XmippWindowUtil.getConstraints(constraints, 0, 3, 1, 1, GridBagConstraints.HORIZONTAL));
			JPanel actionspn = new JPanel(new FlowLayout(FlowLayout.RIGHT));
			actionspn.add(savebt);
			actionspn.add(saveandexitbt);
			add(actionspn, XmippWindowUtil.getConstraints(constraints, 0, 4, 1, 1, GridBagConstraints.HORIZONTAL));

			pack();
			positionx = 0.9f;
			XmippWindowUtil.setLocation(positionx, 0.2f, this);
			setVisible(true);
		}
		catch (Exception e)
		{
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	protected void enableEdition(boolean enable)
	{
		super.enableEdition(enable);

		editfamiliesmi.setEnabled(enable);
		actionsbt.setEnabled(enable);
	}

	public void initMenuBar()
	{
		mb = new JMenuBar();

		// Setting menus

		exportmi = new JMenuItem("Export Particles...", XmippResource.getIcon("export_wiz.gif"));

		exportmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				XmippFileChooser fc = new XmippFileChooser();
				int returnVal = fc.showOpenDialog(TrainingPickerJFrame.this);

				try
				{
					if (returnVal == XmippFileChooser.APPROVE_OPTION)
					{
						File file = fc.getSelectedFile();
						((TrainingPicker) getParticlePicker()).exportParticles(file.getAbsolutePath());
						showMessage("Export successful");
					}
				}
				catch (Exception ex)
				{
					showException(ex);
				}
			}
		});
		filemn.add(importffmi);
		if (ppicker.getFamily().getStep() != FamilyState.Manual)
			importffmi.setEnabled(false);
		filemn.add(exportmi);

		
		
		JMenu windowmn = new JMenu("Window");

		mb.add(filemn);
		mb.add(filtersmn);
		mb.add(windowmn);
		mb.add(helpmn);
		// importffilemi.setText("Import from File...");

		windowmn.add(pmi);
		windowmn.add(ijmi);

		templatesmi = new JMenuItem("Templates");
		editfamiliesmi = new JMenuItem("Edit Families", XmippResource.getIcon("edit.gif"));
		windowmn.add(editfamiliesmi);
		windowmn.add(templatesmi);

		// Setting menu item listeners

		editfamiliesmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				new EditFamiliesJDialog(TrainingPickerJFrame.this, true);

			}
		});

		templatesmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				if (templatesdialog == null)
				{
					templatesdialog = new TemplatesJDialog(TrainingPickerJFrame.this);
					UpdateTemplatesTask.setTemplatesDialog(templatesdialog);
					ParticleToTemplatesTask.setTemplatesDialog(templatesdialog);
				}
				else
				{

					templatesdialog.setVisible(true);
				}

			}
		});
	}



	private void initFamilyPane()
	{

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
		familiescb.setEnabled(ppicker.getMode() == FamilyState.Manual || ppicker.getMode() == FamilyState.ReadOnly);
		if (ppicker.getMode() == FamilyState.Manual && family.getStep() != FamilyState.Manual)
			throw new IllegalArgumentException(
					String.format("Application not enabled for %s mode. Family %s could not be loaded", family.getStep(), family.getName()));

		fieldspn.add(familiescb);

		// Setting color
		initColorPane(family.getColor());
		fieldspn.add(colorpn);

		// Setting slider
		initSizePane();
		fieldspn.add(sizepn);

		centerpickchb = new JCheckBox("Center Particle");
		centerpickchb.setSelected(true);
		fieldspn.add(centerpickchb);

		familypn.add(fieldspn, 0);
		JPanel steppn = new JPanel(new FlowLayout(FlowLayout.LEFT));
		steppn.add(new JLabel("Step:"));
		FamilyState step = family.getStep();
		steplb = new JLabel();
		steppn.add(steplb);

		index = ppicker.getMicrographIndex();

		initThresholdPane();
		steppn.add(thresholdpn);
		actionsbt = XmippWindowUtil.getTextButton("", new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				if (actionsbt.getText().equals(MicrographFamilyState.Autopick.toString()))
					autopick();
				else if (actionsbt.getText().equals(MicrographFamilyState.Correct.toString()))
					correct();
			}
		});
		autopickpercentpn = new JPanel();
		autopickpercentpn.add(new JLabel("Percent to Check"));
		autopickpercenttf = new JFormattedTextField(NumberFormat.getIntegerInstance());
		autopickpercenttf.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				if (autopickpercenttf.getValue() == null)
				{
					JOptionPane.showMessageDialog(TrainingPickerJFrame.this, XmippMessage.getEmptyFieldMsg("Percent to Check"));
					autopickpercenttf.setValue(getFamilyData().getAutopickpercent());
					return;
				}

				int autopickpercent = ((Number) autopickpercenttf.getValue()).intValue();
				getFamilyData().setAutopickpercent(autopickpercent);
				ppicker.setAutopickpercent(autopickpercent);
				ppicker.saveConfig();

			}
		});
		autopickpercenttf.setColumns(3);
		autopickpercentpn.add(autopickpercenttf);

		setStep(step);
		steppn.add(actionsbt);

		steppn.add(autopickpercentpn);

		colorbt.addActionListener(new ColorActionListener());

		familypn.add(steppn, 1);

		familiescb.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				Family family2 = (Family) familiescb.getSelectedItem();
				if (family == family2)
					return;
				// You can only switch between different states for readonly
				// mode. Besides switching will be availabe on manual and
				// readonly modes
				if (family.getStep() != family2.getStep() && ppicker.getMode() != FamilyState.ReadOnly)
				{
					familiescb.setSelectedItem(family);
					JOptionPane.showMessageDialog(TrainingPickerJFrame.this, String
							.format("Application not enabled for %s mode. Family %s could not be loaded", FamilyState.Supervised, family2.getName()));
					return;

				}
				family = family2;
				ppicker.setFamily(family);
				ppicker.saveConfig();

				micrographstb.setRowSelectionInterval(index, index);
				color = (family.getColor());
				colorbt.setIcon(new ColorIcon(color));
				sizesl.setValue(family.getSize());
				updateMicrographsModel();
				micrographstb.getColumnModel().getColumn(micrographsmd.getParticlesPosition()).setHeaderValue(family.getName());
				setThresholdValue(0);
			}
		});

	}

	protected void goToNextStep()
	{

		family.validateNextStep(ppicker);// throws error

		if (family.getStep() == FamilyState.Manual)
			train();

	}

	private void initThresholdPane()
	{
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
		thresholdpn.setVisible(getFamilyData().getState() == MicrographFamilyState.Correct);
		thresholdtf.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				double threshold = Double.parseDouble(thresholdtf.getText()) * 100;
				setThresholdValue(threshold);

			}
		});

		thresholdsl.addChangeListener(new ChangeListener()
		{

			@Override
			public void stateChanged(ChangeEvent e)
			{
				double threshold = (double) thresholdsl.getValue() / 100;
				thresholdtf.setText(String.format("%.2f", threshold));
				setThresholdChanges();
			}
		});
	}

	private void setThresholdValue(double threshold)
	{
		if (Math.abs(threshold) <= 100)
		{
			thresholdsl.setValue((int) threshold);
			setThresholdChanges();
		}
	}

	private void setThresholdChanges()
	{
		// setChanged(true);
		updateMicrographsModel();
		canvas.repaint();
		if (particlesdialog != null)
			loadParticles();

	}

	private void initMicrographsPane()
	{
		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(0, 5, 0, 5);
		constraints.anchor = GridBagConstraints.NORTHWEST;
		micrographpn = new JPanel(new GridBagLayout());
		micrographpn.setBorder(BorderFactory.createTitledBorder("Micrograph"));
		JScrollPane sp = new JScrollPane();
		JPanel ctfpn = new JPanel();
		ctfpn.setBorder(BorderFactory.createTitledBorder(null, "CTF", TitledBorder.CENTER, TitledBorder.BELOW_BOTTOM));
		iconbt = new JButton(Micrograph.getNoImageIcon());
		iconbt.setToolTipText("Load CTF Profile");
		iconbt.setBorderPainted(false);
		iconbt.setContentAreaFilled(false);
		iconbt.setFocusPainted(false);
		iconbt.setOpaque(false);
		iconbt.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				String psd = getMicrograph().getPSD();
				String ctf = getMicrograph().getCTF();
				if (psd != null && ctf != null)
					//ImagesWindowFactory.openCTFWindow(getMicrograph().getPSDImage(), getMicrograph().getCTF(), getMicrograph().getPSD());
					new CTFAnalyzerJFrame(getMicrograph().getPSDImage(), getMicrograph().getCTF(), getMicrograph().getPSD());

			}
		});
		ctfpn.add(iconbt);
		micrographsmd = new MicrographsTableModel(this);
		micrographstb.setModel(micrographsmd);
		formatMicrographsTable();

		sp.setViewportView(micrographstb);
		micrographpn.add(sp, XmippWindowUtil.getConstraints(constraints, 0, 0, 1));
		micrographpn.add(ctfpn, XmippWindowUtil.getConstraints(constraints, 1, 0, 1));
		JPanel infopn = new JPanel();
		manuallb = new JLabel(Integer.toString(ppicker.getManualParticlesNumber(family)));
		autolb = new JLabel(Integer.toString(ppicker.getAutomaticNumber(family, getThreshold())));
		infopn.add(new JLabel("Manual:"));
		infopn.add(manuallb);
		infopn.add(new JLabel("Automatic:"));
		infopn.add(autolb);
		micrographpn.add(infopn, XmippWindowUtil.getConstraints(constraints, 0, 1, 1));
		JPanel buttonspn = new JPanel(new FlowLayout(FlowLayout.LEFT));

		buttonspn.add(resetbt);
		micrographpn.add(buttonspn, XmippWindowUtil.getConstraints(constraints, 0, 2, 2));

	}

	protected void loadMicrograph()
	{
		if (micrographstb.getSelectedRow() == -1)
			return;// Probably from fireTableDataChanged raised
		// is same micrograph??
		if (index == micrographstb.getSelectedRow() && canvas != null && canvas.getIw().isVisible())
			return;
		if (ppicker.isChanged())
			ppicker.saveData(getMicrograph());// Saving changes when switching

		index = micrographstb.getSelectedRow();
		ppicker.getMicrograph().releaseImage();
		ppicker.setMicrograph(ppicker.getMicrographs().get(index));

		ppicker.saveConfig();
		setChanged(false);
		initializeCanvas();
		iconbt.setIcon(ppicker.getMicrograph().getCTFIcon());
		manageAction();
		thresholdpn.setVisible(getFamilyData().getState() == MicrographFamilyState.Correct);
		pack();

		if (particlesdialog != null)
			loadParticles();

	}

	private void manageAction()
	{
		MicrographFamilyData mfd = getFamilyData();
		actionsbt.setText(mfd.getAction());
		boolean isautopick = isAutopick();
		autopickpercentpn.setVisible(isautopick);
		if (isautopick)
		{
			getFamilyData().setAutopickpercent(ppicker.getAutopickpercent());
			autopickpercenttf.setValue(ppicker.getAutopickpercent());
		}
		actionsbt.setVisible(mfd.isActionVisible());
	}

	protected void resetMicrograph()
	{
		ppicker.resetFamilyData(getFamilyData());
		canvas.refreshActive(null);
		updateMicrographsModel();
		setState(MicrographFamilyState.Available);
		
	}

	private void setState(MicrographFamilyState state)
	{
		MicrographFamilyData mfd = getFamilyData();
		mfd.setState(state);
		manageAction();
		ppicker.saveData(getMicrograph());// to keep consistence between files
											// of automatic picker and mines
		setChanged(false);
		thresholdpn.setVisible(state == MicrographFamilyState.Correct);
		updateMicrographsModel();
		pack();
	}

	private boolean isAutopick()
	{
		String action = getFamilyData().getAction();
		if (action == null)
			return false;
		return action.equalsIgnoreCase(MicrographFamilyState.Autopick.toString());

	}

	public MicrographFamilyData getFamilyData()
	{
		return ppicker.getFamilyData();
	}

	private void setStep(FamilyState step)
	{
		MicrographFamilyData mfd = getFamilyData();
		if (micrographsmd != null)
		{
			updateMicrographsModel();
			canvas.repaint();// paints only current class in supervised mode
		}

		steplb.setText(step.toString());
		sizesl.setEnabled(step == FamilyState.Manual);
		sizetf.setEnabled(step == FamilyState.Manual);
		editfamiliesmi.setEnabled(step == FamilyState.Manual);
		manageAction();
		thresholdpn.setVisible(mfd.getState() == MicrographFamilyState.Correct);
		pack();

	}

	protected void initializeCanvas()
	{

		if (canvas == null)
			canvas = new TrainingCanvas(this);
		else
			canvas.updateMicrograph();

		canvas.display();
		updateZoom();
	}

	private void formatMicrographsTable()
	{
		micrographstb.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		micrographstb.getColumnModel().getColumn(0).setPreferredWidth(35);
		micrographstb.getColumnModel().getColumn(1).setPreferredWidth(245);
		micrographstb.getColumnModel().getColumn(2).setPreferredWidth(70);
		micrographstb.getColumnModel().getColumn(3).setPreferredWidth(70);
		micrographstb.setPreferredScrollableViewportSize(new Dimension(420, 304));
		micrographstb.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		if (index != -1)
			micrographstb.setRowSelectionInterval(index, index);
	}

	protected void saveChanges()
	{

		ppicker.saveData();
		setChanged(false);
	}

	void updateFamilyColor()
	{
		color = family.getColor();
		colorbt.setIcon(new ColorIcon(color));
		canvas.repaint();
		ppicker.saveFamilies();
	}

	void updateFamilyComboBox()
	{
		Family item = (Family) familiescb.getSelectedItem();
		DefaultComboBoxModel model = new DefaultComboBoxModel(ppicker.getFamilies().toArray());
		familiescb.setModel(model);
		familiescb.setSelectedItem(item);
		micrographsmd.fireTableStructureChanged();

		formatMicrographsTable();
		pack();
		ppicker.saveFamilies();
	}

	public void addFamily(Family g)
	{
		if (ppicker.existsFamilyName(g.getName()))
			throw new IllegalArgumentException(XmippMessage.getAlreadyExistsGroupNameMsg(g.getName()));
		ppicker.getFamilies().add(g);
		updateFamilyComboBox();
	}

	public void removeFamily(Family family)
	{
		ppicker.removeFamily(family);
		updateFamilyComboBox();
	}

	public void setChanged(boolean changed)
	{
		ppicker.setChanged(changed);
		savemi.setEnabled(changed);
		savebt.setEnabled(changed);
	}

	public void updateMicrographsModel(boolean all)
	{

		

		if (particlesdialog != null)
			loadParticles();

		if (all)
			micrographsmd.fireTableRowsUpdated(0, micrographsmd.getRowCount() - 1);
		else
			micrographsmd.fireTableRowsUpdated(index, index);

		micrographstb.setRowSelectionInterval(index, index);
		manuallb.setText(Integer.toString(ppicker.getManualParticlesNumber(family)));
		autolb.setText(Integer.toString(ppicker.getAutomaticNumber(family, getThreshold())));
		actionsbt.setVisible(getFamilyData().isActionVisible());
	}

	public ParticlePickerCanvas getCanvas()
	{
		return canvas;
	}

	private void train()
	{
		try
		{
			family.goToNextStep(ppicker);// validate and change state if
											// possible
			// setChanged(true);
			setStep(FamilyState.Supervised);// change visual appearance
			ppicker.saveFamilies();

			canvas.setEnabled(false);
			XmippWindowUtil.blockGUI(this, "Training...");

			Thread t = new Thread(new Runnable()
			{

				public void run()
				{
					MicrographFamilyData mfd;
					String args;
					SupervisedParticlePicker sppicker = ((SupervisedParticlePicker) ppicker);
					for (TrainingMicrograph micrograph : ppicker.getMicrographs())
					{
						mfd = micrograph.getFamilyData(family);
						if (!mfd.isEmpty())
						{
							args = sppicker.getBuildInvariantCommandLineArgs(mfd);
							ppicker.runXmippProgram("xmipp_micrograph_automatic_picking", args);
						}
					}
					args = sppicker.getTrainCommandLineArgs();

					ppicker.runXmippProgram("xmipp_micrograph_automatic_picking", args);

					XmippWindowUtil.releaseGUI(getRootPane());
					canvas.setEnabled(true);
				}
			});
			t.start();

		}
		catch (Exception e)
		{
			TrainingPicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	private void autopick()
	{
		setState(MicrographFamilyState.Autopick);

		final String fargs = ((SupervisedParticlePicker) ppicker).getAutopickCommandLineArgs(getFamilyData());
		try
		{
			canvas.setEnabled(false);
			XmippWindowUtil.blockGUI(this, "Autopicking...");
			Thread t = new Thread(new Runnable()
			{

				public void run()
				{
					ppicker.runXmippProgram("xmipp_micrograph_automatic_picking", fargs);
					ppicker.loadAutomaticParticles(getFamilyData());
					
					setState(MicrographFamilyState.Correct);
					canvas.repaint();
					canvas.setEnabled(true);
					XmippWindowUtil.releaseGUI(getRootPane());
				}
			});
			t.start();

		}
		catch (Exception e)
		{
			TrainingPicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	private void correct()
	{
		getFamilyData().deleteBelowThreshold(getThreshold());
		setState(MicrographFamilyState.ReadOnly);
		ppicker.saveAutomaticParticles(getFamilyData());

		try
		{
			canvas.setEnabled(false);
			XmippWindowUtil.blockGUI(this, "Correcting...");

			Thread t = new Thread(new Runnable()
			{
				public void run()
				{
					SupervisedParticlePicker sppicker = ((SupervisedParticlePicker) ppicker);
					String args;
					args = sppicker.getBuildInvariantCommandLineArgs(getFamilyData());
					sppicker.runXmippProgram("xmipp_micrograph_automatic_picking", args);// build
																							// invariants
					args = sppicker.getCorrectCommandLineArgs(getFamilyData());

					ppicker.runXmippProgram("xmipp_micrograph_automatic_picking", args);// correct
					actionsbt.setVisible(false);

					canvas.setEnabled(true);
					XmippWindowUtil.releaseGUI(getRootPane());
					if (index < micrographsmd.getRowCount())
						micrographstb.setRowSelectionInterval(index + 1, index + 1);
				}
			});
			t.start();

		}
		catch (Exception e)
		{
			TrainingPicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public double getThreshold()
	{
		return thresholdsl.getValue() / 100.0;
	}

	@Override
	public List<? extends TrainingParticle> getAvailableParticles()
	{
		return getFamilyData().getAvailableParticles(getThreshold());
	}

	@Override
	public boolean isPickingAvailable(MouseEvent e)
	{
		if (!super.isPickingAvailable(e))
			return false;
		return getFamilyData().isPickingAvailable();
	}

	public String importParticlesFromFile(Format format, String file, float scale, boolean invertx, boolean inverty)
	{
		String result = "";
		if (ppicker.isReviewFile(file))
		{
			result = ppicker.importAllParticles(file, scale, invertx, inverty);
			ppicker.saveData();
		}
		else
			result = importMicrographParticles(format, file, scale, invertx, inverty);
		setChanged(false);
		getCanvas().repaint();
		updateMicrographsModel();
		updateSize(family.getSize());
		canvas.refreshActive(null);
		return result;
	}

	public String importMicrographParticles(Format format, String file, float scale, boolean invertx, boolean inverty)
	{
		family.initTemplates();
		String filename = Micrograph.getName(file, 1);
		// validating you want use this file for this micrograph with different
		// name
		if (!filename.equals(getMicrograph().getName()))
		{
			String msg = String.format("Are you sure you want to import data from file\n%s to micrograph %s ?", file, getMicrograph().getName());
			int result = JOptionPane.showConfirmDialog(this, msg);
			if (result != JOptionPane.YES_OPTION)
				return null;
		}
		MicrographFamilyData mfd = getFamilyData();
		mfd.reset();
		String result = ((ManualParticlePicker) ppicker).importParticlesFromFile(file, format, mfd.getMicrograph(), scale, invertx, inverty);
		ppicker.saveData(getMicrograph());
		return result;
	}

	@Override
	protected void openHelpURl()
	{
		XmippWindowUtil.openURI("http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Micrograph_particle_picking_v3");

	}



	public void updateSize(int size)
	{
		try
		{
			ppicker.resetParticleImages();
			super.updateSize(size);
			ppicker.updateTemplates();
			
		}
		catch (Exception e)
		{
			String msg = (e.getMessage() != null) ? e.getMessage() : XmippMessage.getUnexpectedErrorMsg();
			XmippDialog.showError(this, msg);
		}
	}

	public boolean isCenterParticle()
	{
		return centerpickchb.isSelected();
	}


	@Override
	public String importParticles(Format format, String dir, float scale, boolean invertx, boolean inverty)
	{
		String result = "";

		if (new File(dir).isDirectory())
		{
			((ManualParticlePicker) ppicker).importParticlesFromFolder(dir, format, scale, invertx, inverty);
			getCanvas().repaint();
			updateMicrographsModel(true);
			getCanvas().refreshActive(null);
		}
		else
			// only can choose file if TrainingPickerJFrame instance
			result = importParticlesFromFile(format, dir, scale, invertx, inverty);
		
		return result;

	}



	@Override
	public ParticlesJDialog initParticlesJDialog()
	{
		return new ParticlesJDialog(this);
	}

	
}
