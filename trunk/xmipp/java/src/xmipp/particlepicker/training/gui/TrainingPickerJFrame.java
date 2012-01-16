package xmipp.particlepicker.training.gui;

import ij.IJ;
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
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import xmipp.particlepicker.Constants;
import xmipp.particlepicker.Family;
import xmipp.particlepicker.ParticlePickerCanvas;
import xmipp.particlepicker.ParticlePickerJFrame;
import xmipp.particlepicker.WindowUtils;
import xmipp.particlepicker.training.model.FamilyState;
import xmipp.particlepicker.training.model.MicrographFamilyData;
import xmipp.particlepicker.training.model.MicrographFamilyState;
import xmipp.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.particlepicker.training.model.TrainingMicrograph;
import xmipp.particlepicker.training.model.TrainingParticle;
import xmipp.particlepicker.training.model.TrainingPicker;
import xmipp.jni.Program;


public class TrainingPickerJFrame extends ParticlePickerJFrame 
{

	private TrainingCanvas canvas;
	private JMenuBar mb;
	private JComboBox familiescb;
	private TrainingPicker ppicker;
	private Color color;
	private JPanel familypn;
	private JPanel micrographpn;
	private JTable micrographstb;

	private MicrographsTableModel micrographsmd;
	private TrainingMicrograph micrograph;
	private JButton nextbt;
	private JButton colorbt;
	private double positionx;
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
	
	private JMenuItem exportmi;




	@Override
	public TrainingPicker getParticlePicker()
	{
		return ppicker;
	}

	public Family getFamily()
	{
		return family;
	}

	public TrainingPickerJFrame(TrainingPicker picker)
	{
		super(picker);
		this.ppicker = picker;
		initComponents();
	}

	public TrainingMicrograph getMicrograph()
	{
		return micrograph;
	}

	private void initComponents()
	{
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
		add(familypn, WindowUtils.getConstraints(constraints, 0, 1, 3));

		initSymbolPane();
		add(symbolpn, WindowUtils.getConstraints(constraints, 0, 2, 3));

		initMicrographsPane();
		add(micrographpn, WindowUtils.getConstraints(constraints, 0, 3, 3));

		pack();
		positionx = 0.995;
		WindowUtils.centerScreen(positionx, 0.25, this);
		setVisible(true);
	}

	public void initMenuBar()
	{
		mb = new JMenuBar();

		// Setting menus
		JMenu filemn = new JMenu("File");
		JMenu windowmn = new JMenu("Window");
		JMenu helpmn = new JMenu("Help");
		mb.add(filemn);
		mb.add(filtersmn);
		mb.add(windowmn);
		mb.add(helpmn);

		// Setting menu items
		savemi.setEnabled(ppicker.isChanged());
		filemn.add(savemi);

		exportmi = new JMenuItem("Export Particles");
		filemn.add(exportmi);
		windowmn.add(pmi);
		windowmn.add(ijmi);
		editfamiliesmn = new JMenuItem("Edit Families");
		windowmn.add(editfamiliesmn);

		helpmn.add(hcontentsmi);

		// Setting menu item listeners

	

		exportmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				ppicker.exportData(family);
				JOptionPane.showMessageDialog(TrainingPickerJFrame.this, "Export successful");
			}
		});

		

		editfamiliesmn.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				new EditFamiliesJDialog(TrainingPickerJFrame.this, true);

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
		familiescb.setEnabled(ppicker.getMode() != FamilyState.Review);

		family = (Family) familiescb.getSelectedItem();
		if (ppicker.getMode() == FamilyState.Manual && family.getStep() != FamilyState.Manual)
			throw new IllegalArgumentException(String.format("Application not enabled for %s mode. Family %s could not be loaded", family.getStep(), family.getName()));

		fieldspn.add(familiescb);

		// Setting color
		color = family.getColor();
		fieldspn.add(new JLabel("Color:"));
		colorbt = new JButton();
		colorbt.setIcon(new ColorIcon(color));
		colorbt.setBorderPainted(false);
		fieldspn.add(colorbt);

		// Setting slider
		initSizePane();
		fieldspn.add(sizepn);

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
		nextbt.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				try
				{
					family.validateNextStep(ppicker);
				}
				catch (Exception ex)
				{
					JOptionPane.showMessageDialog(TrainingPickerJFrame.this, ex.getMessage());
					return;
				}
				if (family.getStep() == FamilyState.Manual)
				{
					int result = JOptionPane.showConfirmDialog(TrainingPickerJFrame.this, "Data provided from manual picking will be used to train software and start Supervised mode." + "\nManual mode will be disabled. Are you sure you want to continue?", "Message", JOptionPane.YES_NO_OPTION);

					if (result == JOptionPane.NO_OPTION)
						return;
					train();

				}
				else if (family.getStep() == FamilyState.Supervised)
				{
					int result = JOptionPane.showConfirmDialog(TrainingPickerJFrame.this, "Model generated during Supervised mode will be dismissed." + "\nAre you sure you want to continue?", "Message", JOptionPane.YES_NO_OPTION);

					if (result == JOptionPane.NO_OPTION)
						return;
					family.goToPreviousStep();
					((SupervisedParticlePicker) ppicker).resetModel(family);
					setStep(FamilyState.Manual);

				}
			}
		});
		familypn.add(steppn, 1);
		
		familiescb.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				Family family2 = (Family) familiescb.getSelectedItem();
				if (ppicker.getMode() == FamilyState.Manual && family2.getStep() == FamilyState.Supervised)
				{
					familiescb.setSelectedItem(family);
					JOptionPane.showMessageDialog(TrainingPickerJFrame.this, String.format("Application not enabled for %s mode. Family %s could not be loaded", FamilyState.Supervised, family2.getName()));
					return;
				}
				if (family.getStep() != family2.getStep())
				{
					int result = JOptionPane.showConfirmDialog(TrainingPickerJFrame.this, String.format("Selecting family %s will take you to %s mode." + "\nAre you sure you want to continue?", family2.getName(), family2.getStep().toString()), "Message", JOptionPane.YES_NO_OPTION);
					if (result == JOptionPane.NO_OPTION)
					{
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
				micrographstb.getColumnModel().getColumn(micrographsmd.getParticlesPosition()).setHeaderValue(family.getName());
			}
		});
		actionsbt.addActionListener(new ActionListener()
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

	}

	class ColorActionListener implements ActionListener
	{
		JColorChooser colorChooser;

		@Override
		public void actionPerformed(ActionEvent e)
		{
			// Set up the dialog that the button brings up.
			colorChooser = new JColorChooser();
			JDialog dialog = JColorChooser.createDialog(colorbt, "Pick a Color", true, // modal
					colorChooser, new ActionListener()
					{

						@Override
						public void actionPerformed(ActionEvent e)
						{
							family.setColor(ColorActionListener.this.colorChooser.getColor());
							updateFamilyColor();
						}
					}, // OK button handler
					null); // no CANCEL button handler
			WindowUtils.centerScreen(positionx, 0.25, dialog);
			dialog.setVisible(true);
		}
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
				if (Math.abs(threshold) <= 100)
				{
					thresholdsl.setValue((int) threshold);
					setThresholdChanges();
				}

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

	private void setThresholdChanges()
	{
		updateMicrographsModel();
		canvas.repaint();
		actionsbt.setVisible(getFamilyData().isActionAvailable(getThreshold()));
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
		ctfpn.setBorder(BorderFactory.createTitledBorder(null, "CTF", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.BELOW_BOTTOM));
		iconlb = new JLabel();
		ctfpn.add(iconlb);
		micrographsmd = new MicrographsTableModel(this);
		micrographstb = new JTable(micrographsmd);
		micrographstb.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		micrographstb.getColumnModel().getColumn(0).setPreferredWidth(35);
		micrographstb.getColumnModel().getColumn(1).setPreferredWidth(245);
		micrographstb.getColumnModel().getColumn(2).setPreferredWidth(70);
		micrographstb.getColumnModel().getColumn(3).setPreferredWidth(70);
		micrographstb.setPreferredScrollableViewportSize(new Dimension(420, 304));
		micrographstb.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

		sp.setViewportView(micrographstb);
		micrographpn.add(sp, WindowUtils.getConstraints(constraints, 0, 0, 1));
		micrographpn.add(ctfpn, WindowUtils.getConstraints(constraints, 1, 0, 1));
		JPanel infopn = new JPanel();
		manuallb = new JLabel(Integer.toString(ppicker.getManualParticlesNumber(family)));
		autolb = new JLabel(Integer.toString(ppicker.getAutomaticNumber(family, getThreshold())));
		infopn.add(new JLabel("Manual:"));
		infopn.add(manuallb);
		infopn.add(new JLabel("Automatic:"));
		infopn.add(autolb);
		micrographpn.add(infopn, WindowUtils.getConstraints(constraints, 0, 1, 1));
		JPanel buttonspn = new JPanel(new FlowLayout(FlowLayout.LEFT));
		resetbt = new JButton("Reset");
		buttonspn.add(resetbt);
		micrographpn.add(buttonspn, WindowUtils.getConstraints(constraints, 0, 2, 2));
		resetbt.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				ppicker.resetFamilyData(getFamilyData());
				setState(MicrographFamilyState.Available);
				canvas.setActive(null);
				updateMicrographsModel();
				setChanged(true);
			}
		});
		micrographstb.getSelectionModel().addListSelectionListener(new ListSelectionListener()
		{

			@Override
			public void valueChanged(ListSelectionEvent e)
			{
				if (e.getValueIsAdjusting())
					return;
				if (TrainingPickerJFrame.this.micrographstb.getSelectedRow() == -1)
					return;// Probably from fireTableDataChanged raised
				index = TrainingPickerJFrame.this.micrographstb.getSelectedRow();
				// by me.
				micrograph.releaseImage();
				micrograph = (TrainingMicrograph) ppicker.getMicrographs().get(index);

				initializeCanvas();
				TrainingPickerJFrame.this.iconlb.setIcon(micrograph.getCTFIcon());
				actionsbt.setText(getFamilyData().getAction());
				actionsbt.setVisible(getFamilyData().isActionAvailable(getThreshold()));
				thresholdpn.setVisible(getFamilyData().getState() == MicrographFamilyState.Correct);
				pack();
				saveChanges();// Saving changes when switching micrographs, by
								// Coss suggestion
				if (particlesdialog != null)
					loadParticles();
			}
		});
		micrographstb.getSelectionModel().setSelectionInterval(index, index);

	}

	private void setState(MicrographFamilyState state)
	{
		getFamilyData().setState(state);
		actionsbt.setText(getFamilyData().getAction());
		saveChanges();// to keep consistence between files of automatic picker
						// and mines
		thresholdpn.setVisible(state == MicrographFamilyState.Correct);
		updateMicrographsModel();
		pack();
	}

	MicrographFamilyData getFamilyData()
	{
		return micrograph.getFamilyData(family);
	}

	private void setStep(FamilyState step)
	{
		if (micrographsmd != null)
		{
			updateMicrographsModel();
			canvas.repaint();// paints only current class in supervised mode
		}
		if (step == FamilyState.Manual)
			nextbt.setText("Go To " + TrainingPicker.nextStep(step).toString());
		else if (step == FamilyState.Supervised)
			nextbt.setText("Go To " + TrainingPicker.previousStep(step).toString());
		else
			nextbt.setVisible(false);
		steplb.setText(step.toString());
		sizesl.setEnabled(step == FamilyState.Manual);
		sizetf.setEnabled(step == FamilyState.Manual);
		editfamiliesmn.setEnabled(step == FamilyState.Manual);
		actionsbt.setText(getFamilyData().getAction());
		actionsbt.setVisible(getFamilyData().isActionAvailable(getThreshold()));
		thresholdpn.setVisible(getFamilyData().getState() == MicrographFamilyState.Correct);
		pack();
	}

	protected void initializeCanvas()
	{
		if (canvas == null)
		{
			canvas = new TrainingCanvas(this);
			ImageWindow iw = new ImageWindow(micrograph.getImagePlus(), canvas);
			iw.setTitle(micrograph.getName());
			if(!ppicker.getFilters().isEmpty())
				IJ.runMacro(ppicker.getFiltersMacro());
		}
		else
			canvas.updateMicrograph();
	}

	protected void saveChanges()
	{
		ppicker.saveData();
		setChanged(false);
	}

	void updateFamilyColor()
	{
		setChanged(true);
		color = family.getColor();
		colorbt.setIcon(new ColorIcon(color));
		canvas.repaint();
	}

	void updateFamilyComboBox()
	{
		setChanged(true);
		Family item = (Family) familiescb.getSelectedItem();
		DefaultComboBoxModel model = new DefaultComboBoxModel(ppicker.getFamilies().toArray());
		familiescb.setModel(model);
		familiescb.setSelectedItem(item);
		pack();
	}



	public void addFamily(Family g)
	{
		if (ppicker.existsFamilyName(g.getName()))
			throw new IllegalArgumentException(Constants.getAlreadyExistsGroupNameMsg(g.getName()));
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
	}

	public void updateMicrographsModel()
	{
		super.updateMicrographsModel();
		micrographsmd.fireTableRowsUpdated(index, index);
		micrographstb.setRowSelectionInterval(index, index);
		manuallb.setText(Integer.toString(ppicker.getManualParticlesNumber(family)));
		autolb.setText(Integer.toString(ppicker.getAutomaticNumber(family, getThreshold())));
		actionsbt.setVisible(getFamilyData().isActionAvailable(getThreshold()));
	}

	public ParticlePickerCanvas getCanvas()
	{
		return canvas;
	}

	private void train()
	{

		
		family.goToNextStep(ppicker);//validate and change state if posible
		//setChanged(true);
		setStep(FamilyState.Supervised);//change visual appearance
		saveChanges();//persist changes
		String args;
		for (TrainingMicrograph micrograph : ppicker.getMicrographs())
		{
			if (!micrograph.getFamilyData(family).isEmpty())
			{
				
				args = String.format("-i %s --particleSize %s --model %s --outputRoot %s --mode train %s", micrograph.getFile(),// -i
						family.getSize(), // --particleSize
						ppicker.getOutputPath(family.getName()),// --model
						ppicker.getOutputPath(micrograph.getName()), // --outputRoot
						family.getName() + "@" + ppicker.getOutputPath(micrograph.getOFilename()));// train
				// parameter
				if (((SupervisedParticlePicker) ppicker).isFastMode())
					args += " --fast";
				if (((SupervisedParticlePicker) ppicker).isIncore())
					args += " --in_core";
				final String fargs = args;
				try
				{
					final InfiniteProgressPanel glassPane = new InfiniteProgressPanel("training picker...");
					final Component previousGlassPane = TrainingPickerJFrame.this.getRootPane().getGlassPane();
					canvas.setEnabled(false);
					TrainingPickerJFrame.this.getRootPane().setGlassPane(glassPane);
					glassPane.start();

					Thread t = new Thread(new Runnable()
					{

						public void run()
						{

							try
							{

								Program.runByName("xmipp_micrograph_automatic_picking", fargs);
							}
							catch (Exception e)
							{
								TrainingPicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
								throw new IllegalArgumentException(e);
							}
							glassPane.stop();
							TrainingPickerJFrame.this.getRootPane().setGlassPane(previousGlassPane);
							canvas.setEnabled(true);
							int next = ppicker.getNextFreeMicrograph(family);
							if (next != -1)
								micrographstb.setRowSelectionInterval(next, next);

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
		}
	}

	private void autopick()
	{
		setState(MicrographFamilyState.Autopick);
		String args;
		args = String.format("-i %s --particleSize %s --model %s --outputRoot %s --mode try --thr %s", micrograph.getFile(),// -i
				family.getSize(), // --particleSize
				ppicker.getOutputPath(family.getName()),// --model
				ppicker.getOutputPath(micrograph.getName()),// --outputRoot
				((SupervisedParticlePicker) ppicker).getThreads()// --thr
		);

		if (((SupervisedParticlePicker) ppicker).isFastMode())
			args += " --fast";
		if (((SupervisedParticlePicker) ppicker).isIncore())
			args += " --in_core";
		final String fargs = args;
		try
		{
			final InfiniteProgressPanel glassPane = new InfiniteProgressPanel("autopicking...");
			final Component previousGlassPane = TrainingPickerJFrame.this.getRootPane().getGlassPane();
			canvas.setEnabled(false);
			TrainingPickerJFrame.this.getRootPane().setGlassPane(glassPane);
			glassPane.start();

			Thread t = new Thread(new Runnable()
			{

				public void run()
				{
					try
					{

						Program.runByName("xmipp_micrograph_automatic_picking", fargs);
					}
					catch (Exception e)
					{
						TrainingPicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
						throw new IllegalArgumentException(e);
					}
					glassPane.stop();
					TrainingPickerJFrame.this.getRootPane().setGlassPane(previousGlassPane);
					ppicker.loadAutomaticParticles(micrograph);
					setState(MicrographFamilyState.Correct);
					canvas.repaint();
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

	private void correct()
	{
		getFamilyData().deleteBelowThreshold(getThreshold());
		setState(MicrographFamilyState.ReadOnly);
		ppicker.persistAutomaticParticles(getFamilyData());

		String args = String.format("-i %s --particleSize %s --model %s --outputRoot %s --mode train ", micrograph.getFile(),// -i
				family.getSize(), // --particleSize
				ppicker.getOutputPath(family.getName()),// --model
				ppicker.getOutputPath(micrograph.getName())// --outputRoot
		);

		if (micrograph.getFamilyData(family).getManualParticles().size() > 0)
			args += family.getName() + "@" + ppicker.getOutputPath(micrograph.getOFilename());
		if (((SupervisedParticlePicker) ppicker).isFastMode())
			args += " --fast";
		if (((SupervisedParticlePicker) ppicker).isIncore())
			args += " --in_core";
		final String fargs = args;
		try
		{
			final InfiniteProgressPanel glassPane = new InfiniteProgressPanel("correcting...");
			final Component previousGlassPane = TrainingPickerJFrame.this.getRootPane().getGlassPane();
			canvas.setEnabled(false);
			TrainingPickerJFrame.this.getRootPane().setGlassPane(glassPane);
			glassPane.start();
			Thread t = new Thread(new Runnable()
			{
				public void run()
				{
					try
					{
						Program.runByName("xmipp_micrograph_automatic_picking", fargs);
					}
					catch (Exception e)
					{
						TrainingPicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
						throw new IllegalArgumentException(e);
					}

					glassPane.stop();

					TrainingPickerJFrame.this.getRootPane().setGlassPane(previousGlassPane);
					int next = ppicker.getNextFreeMicrograph(family);
					if (next != -1)
						micrographstb.setRowSelectionInterval(next, next);
					else
						actionsbt.setVisible(false);
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

	public double getThreshold()
	{
		return thresholdsl.getValue() / 100.0;
	}

	@Override
	public List<? extends TrainingParticle> getParticles()
	{
		return getFamilyData().getParticles();
	}

	@Override
	public boolean isPickingAvailable(MouseEvent e)
	{
	    if(!super.isPickingAvailable(e))
	    	return false;
		return getFamilyData().isPickingAvailable();
	}

	@Override
	public void changeShapes()
	{
		canvas.repaint();
		
	}

}
