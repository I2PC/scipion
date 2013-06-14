package xmipp.viewer.particlepicker.training.gui;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.text.NumberFormat;
import java.util.List;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
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
import xmipp.viewer.particlepicker.Format;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.ParticlesJDialog;
import xmipp.viewer.particlepicker.SingleParticlePicker;
import xmipp.viewer.particlepicker.training.model.MicrographState;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.TrainingMicrograph;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;

public class SingleParticlePickerJFrame extends ParticlePickerJFrame
{

	private SingleParticlePickerCanvas canvas;
	private JMenuBar mb;
	private SingleParticlePicker ppicker;
	private JPanel micrographpn;
	private MicrographsTableModel micrographsmd;

	private float positionx;
	private JButton iconbt;
	private JLabel manuallb;
	private JLabel autolb;
	private JSlider thresholdsl;
	private JPanel thresholdpn;
	private JFormattedTextField thresholdtf;

	private JMenuItem advancedoptionsmi;
	AdvancedOptionsJDialog optionsdialog;
	private JCheckBox centerparticlebt;
	private JMenuItem exportmi;
	private JCheckBox autopickchb;
	private JPanel sppickerpn;
	private JLabel autopicklb;

	@Override
	public SingleParticlePicker getParticlePicker()
	{
		return ppicker;
	}

	public SingleParticlePickerJFrame(SingleParticlePicker picker)
	{

		super(picker);
		try
		{
			this.ppicker = picker;
			initComponents();
			
			
		}
		catch (IllegalArgumentException ex)
		{
			close();
			throw ex;
		}
	}

	public boolean isCenterParticle()
	{
		return centerparticlebt.isSelected();
	}

	@Override
	public ParticlesJDialog initParticlesJDialog()
	{
		return new ParticlesJDialog(this);
	}

	public TrainingMicrograph getMicrograph()
	{
		return ppicker.getMicrograph();
	}

	@Override
	protected void openHelpURl()
	{
		XmippWindowUtil.openURI("http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Micrograph_particle_picking_v3");

	}

	public double getThreshold()
	{
		if(thresholdsl == null)
			return 0;
		return thresholdsl.getValue() / 100.0;
	}

	@Override
	public List<? extends TrainingParticle> getAvailableParticles()
	{
		return getMicrograph().getAvailableParticles(getThreshold());
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
		updateSize(ppicker.getSize());
		canvas.refreshActive(null);
		return result;
	}

	public String importMicrographParticles(Format format, String file, float scale, boolean invertx, boolean inverty)
	{

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
		ppicker.initTemplates();
		getMicrograph().reset();
		String result = ppicker.importParticlesFromFile(file, format, getMicrograph(), scale, invertx, inverty);
		ppicker.saveData(getMicrograph());
		return result;
	}

	public void updateSize(int size)
	{
		try
		{
			ppicker.resetParticleImages();
			super.updateSize(size);

		}
		catch (Exception e)
		{
			String msg = (e.getMessage() != null) ? e.getMessage() : XmippMessage.getUnexpectedErrorMsg();
			XmippDialog.showError(this, msg);
		}
	}

	@Override
	public String importParticles(Format format, String dir, float scale, boolean invertx, boolean inverty)
	{
		String result = "";

		if (new File(dir).isDirectory())
		{
			ppicker.importParticlesFromFolder(dir, format, scale, invertx, inverty);
			getCanvas().repaint();
			updateMicrographsModel(true);
			getCanvas().refreshActive(null);
		}
		else
			// only can choose file if TrainingPickerJFrame instance
			result = importParticlesFromFile(format, dir, scale, invertx, inverty);
		return result;

	}

	protected void initializeCanvas()
	{

		if (canvas == null)
			canvas = new SingleParticlePickerCanvas(this);
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
		int index = ppicker.getMicrographIndex();
		if (index != -1)
			micrographstb.setRowSelectionInterval(index, index);
	}

	void updateFamilyColor()
	{
		color = ppicker.getColor();
		colorbt.setIcon(new ColorIcon(color));
		canvas.repaint();
		ppicker.saveConfig();
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

		int index = ppicker.getMicrographIndex();
		if (all)
			micrographsmd.fireTableRowsUpdated(0, micrographsmd.getRowCount() - 1);
		else
			micrographsmd.fireTableRowsUpdated(index, index);

		micrographstb.setRowSelectionInterval(index, index);
		manuallb.setText(Integer.toString(ppicker.getManualParticlesNumber()));
		autolb.setText(Integer.toString(ppicker.getAutomaticParticlesNumber(getThreshold())));
	}

	public ParticlePickerCanvas getCanvas()
	{
		return canvas;
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

			initToolBar();
			centerparticlebt = new JCheckBox("Center Particle", true);
			tb.add(centerparticlebt);
			add(tb, XmippWindowUtil.getConstraints(constraints, 0, 0, 2, 1, GridBagConstraints.HORIZONTAL));

			add(new JLabel("Shape:"), XmippWindowUtil.getConstraints(constraints, 0, 1));
			initShapePane();
			add(shapepn, XmippWindowUtil.getConstraints(constraints, 1, 1));

			if (ppicker.getMode() != Mode.Review)
			{
				autopicklb = new JLabel("Autopick:");
				constraints.insets = new Insets(30, 5, 30, 5);
				add(autopicklb, XmippWindowUtil.getConstraints(constraints, 0, 3));
				constraints.insets = new Insets(0, 5, 0, 5);
				initSupervisedPickerPane();
				add(sppickerpn, XmippWindowUtil.getConstraints(constraints, 1, 3, 1, 1, GridBagConstraints.HORIZONTAL));
				enableSupervised(ppicker.getMode() == Mode.Supervised);
			}
			initMicrographsPane();
			add(micrographpn, XmippWindowUtil.getConstraints(constraints, 0, 4, 2, 1, GridBagConstraints.HORIZONTAL));
			JPanel actionspn = new JPanel(new FlowLayout(FlowLayout.RIGHT));

			actionspn.add(savebt);
			actionspn.add(saveandexitbt);
			add(actionspn, XmippWindowUtil.getConstraints(constraints, 0, 5, 2, 1, GridBagConstraints.HORIZONTAL));
			if (ppicker.getMode() == Mode.ReadOnly)
				enableEdition(false);
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

	private void initSupervisedPickerPane()
	{
		sppickerpn = new JPanel(new FlowLayout(FlowLayout.LEFT));

		// sppickerpn.setBorder(BorderFactory.createTitledBorder("Supervised Picker"));

		JPanel steppn = new JPanel(new FlowLayout(FlowLayout.LEFT));
		autopickchb = new JCheckBox();
		autopickchb.setSelected(ppicker.getMode() == Mode.Supervised);

		autopickchb.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				try
				{
					boolean isautopick = autopickchb.isSelected();

					if (isautopick)
					{
						ppicker.setMode(Mode.Supervised);
						ppicker.train(SingleParticlePickerJFrame.this);

					}
					else
					{
						ppicker.setMode(Mode.Manual);
						updateMicrographsModel();
						canvas.repaint();
					}
					ppicker.saveConfig();
					enableSupervised(isautopick);
					
				}
				catch (Exception ex)
				{
					ex.printStackTrace();
					XmippDialog.showError(SingleParticlePickerJFrame.this, ex.getMessage());
					autopickchb.setSelected(false);
					return;
				}
			}
		});
		sppickerpn.add(autopickchb);
		if (ppicker.getMode() == Mode.Supervised)
			sppickerpn.add(steppn);
		initThresholdPane();
		sppickerpn.add(thresholdpn);

	}

	protected void enableSupervised(boolean selected)
	{

		thresholdpn.setVisible(selected);
		sizesl.setEnabled(!selected);//not really, if there is some micrograph in sup mode size cannot be changed
		sizetf.setEnabled(!selected);
		importmi.setEnabled(!selected);
		if(optionsdialog != null && optionsdialog.isVisible())
			optionsdialog.enableOptions();
		pack();

	}

	public boolean isSupervised()
	{
		return autopickchb.isSelected();
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
				int returnVal = fc.showOpenDialog(SingleParticlePickerJFrame.this);

				try
				{
					if (returnVal == XmippFileChooser.APPROVE_OPTION)
					{
						File file = fc.getSelectedFile();
						ppicker.exportParticles(file.getAbsolutePath());
						showMessage("Export successful");
					}
				}
				catch (Exception ex)
				{
					showException(ex);
				}
			}
		});
		filemn.add(importmi);
		if (ppicker.getMode() != Mode.Manual)
			importmi.setEnabled(false);
		filemn.add(exportmi);
		filemn.add(exitmi);

		JMenu windowmn = new JMenu("Window");

		mb.add(filemn);
		mb.add(filtersmn);
		mb.add(windowmn);
		mb.add(helpmn);
		// importffilemi.setText("Import from File...");

		windowmn.add(pmi);
		windowmn.add(ijmi);

		advancedoptionsmi = new JMenuItem("Advanced Options");
		windowmn.add(advancedoptionsmi);

		// Setting menu item listeners
		advancedoptionsmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				loadAdvancedOptions();

			}
		});
	}

	public void loadAdvancedOptions()
	{
		try
		{
			if (optionsdialog == null)
				optionsdialog = new AdvancedOptionsJDialog(SingleParticlePickerJFrame.this);
			else
			{
				optionsdialog.setVisible(true);
				optionsdialog.enableOptions();
			}
		}
		catch (Exception e)
		{
			XmippDialog.showError(this, e.getMessage());
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
		thresholdtf.setValue(0);
		thresholdpn.add(thresholdtf);
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
		micrographpn.setBorder(BorderFactory.createTitledBorder("Micrographs"));
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
		manuallb = new JLabel(Integer.toString(ppicker.getManualParticlesNumber()));
		autolb = new JLabel(Integer.toString(ppicker.getAutomaticParticlesNumber(getThreshold())));
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
		int index = ppicker.getMicrographIndex();
		if (index == micrographstb.getSelectedRow() && canvas != null && canvas.getIw().isVisible())
			return;
		if (ppicker.isChanged())
			ppicker.saveData(getMicrograph());// Saving changes when switching
		index = micrographstb.getSelectedRow();
		TrainingMicrograph next = ppicker.getMicrographs().get(index);
		
		tryCorrectAndAutopick(getMicrograph(), next);

		ppicker.getMicrograph().releaseImage();
		ppicker.setMicrograph(next);
		initializeCanvas();
		if (ppicker.getMode() == Mode.Supervised)
			resetbt.setEnabled(getMicrograph().getState() != MicrographState.Manual);
		setChanged(false);
		
		iconbt.setIcon(ppicker.getMicrograph().getCTFIcon());

		pack();

		if (particlesdialog != null)
			loadParticles();

	}

	private void tryCorrectAndAutopick(TrainingMicrograph current, TrainingMicrograph next)
	{
		if(current.equals(next))//app just started
			return;
		boolean isautopick = ppicker.getMode() == Mode.Supervised && next.getState() == MicrographState.Available;
		if (ppicker.getMode() == Mode.Supervised && current.getState() == MicrographState.Supervised)
		{
			boolean iscorrect = XmippDialog.showQuestion(this, "Would you like to correct training with added and deleted particles?");
			if (iscorrect)
			{
				ppicker.correctAndAutopick(this, current, next);
				isautopick = false;
			}

		}
		if (isautopick)//if not done before
			ppicker.autopick(this, next);
	}

	protected void resetMicrograph()
	{
		getMicrograph().reset();
		ppicker.initUpdateTemplates();
		canvas.refreshActive(null);
		updateMicrographsModel();
		setState(MicrographState.Available);
	}

	private void setState(MicrographState state)
	{
		getMicrograph().setState(state);
		ppicker.saveData(getMicrograph());// to keep consistence between files
											// of automatic picker and mines
		setChanged(false);
		updateMicrographsModel();
		pack();
	}
	
	public boolean isPickingAvailable(MouseEvent e)
	{
		if(!super.isPickingAvailable())
			return false;
		return getMicrograph().getState() != MicrographState.Corrected;
	}

}
