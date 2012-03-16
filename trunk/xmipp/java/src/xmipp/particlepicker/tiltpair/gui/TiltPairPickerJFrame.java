package xmipp.particlepicker.tiltpair.gui;

import ij.IJ;
import ij.WindowManager;
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
import java.io.File;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JColorChooser;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import xmipp.particlepicker.Family;
import xmipp.particlepicker.ParticlePickerCanvas;
import xmipp.particlepicker.ParticlePickerJFrame;
import xmipp.particlepicker.tiltpair.model.TiltPairPicker;
import xmipp.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.particlepicker.training.gui.ColorIcon;
import xmipp.particlepicker.training.model.TrainingParticle;
import xmipp.utils.XmippIJUtil;




public class TiltPairPickerJFrame extends ParticlePickerJFrame 
{

	private TiltPairPicker pppicker;
	private JMenuBar mb;
	private Color color;
	private JPanel particlespn;
	private JPanel micrographpn;
	private JTable micrographstb;
	private UntiltedMicrographCanvas canvas;
	private MicrographPairsTableModel micrographsmd;
	private UntiltedMicrograph untiltedmic;
	private JButton colorbt;
	private float position;
	private int index;
	private JButton resetbt;
	private JLabel upslb;
	private TiltedMicrographCanvas tiltedcanvas;
	private JCheckBoxMenuItem anglesmi;

	public TiltPairPicker getParticlePicker()
	{
		return pppicker;
	}
	


	public TiltPairPickerJFrame(TiltPairPicker picker)
	{
		super(picker);
		pppicker = picker;
		initComponents();
	}

	private void initComponents()
	{
		setResizable(false);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Xmipp Particle Pair Picker");
		initMenuBar();
		setJMenuBar(mb);

		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(0, 5, 0, 5);
		constraints.anchor = GridBagConstraints.WEST;
		setLayout(new GridBagLayout());

		initParticlesPane();
		add(particlespn, XmippIJUtil.getConstraints(constraints, 0, 1, 3));
		initMicrographsPane();
		add(micrographpn, XmippIJUtil.getConstraints(constraints, 0, 2, 3));

		pack();
		position = 0.95f;
		XmippIJUtil.setLocation(position, 0.5f, this);
		setVisible(true);
	}
	


	public void initMenuBar()
	{
		mb = new JMenuBar();

		// Setting menus
		JMenu viewmn = new JMenu("View");
		JMenu helpmn = new JMenu("Help");
		mb.add(filemn);
		mb.add(filtersmn);
		mb.add(viewmn);
		mb.add(helpmn);

		
		
		anglesmi = new JCheckBoxMenuItem("Angles");
		anglesmi.setEnabled(false);
		anglesmi.addActionListener(new ActionListener()
		{
			
			@Override
			public void actionPerformed(ActionEvent e)
			{
				canvas.repaint();
				tiltedcanvas.repaint();
			}
		});
		viewmn.add(anglesmi);
		viewmn.add(pmi);
		viewmn.add(ijmi);
		helpmn.add(hcontentsmi);
	}




	private void initParticlesPane()
	{
		particlespn = new JPanel();
		GridLayout gl = new GridLayout(2, 1);
		particlespn.setLayout(gl);

		particlespn.setBorder(BorderFactory.createTitledBorder("Particles"));

		JPanel fieldspn = new JPanel(new FlowLayout(FlowLayout.LEFT));

		// Setting color
		color = pppicker.getColor();
		fieldspn.add(new JLabel("Color:"));
		colorbt = new JButton();
		colorbt.setIcon(new ColorIcon(color));
		colorbt.setBorderPainted(false);
		fieldspn.add(colorbt);

		// Setting slider
		initSizePane();
		fieldspn.add(sizepn);

		particlespn.add(fieldspn, 0);
		initSymbolPane();
		particlespn.add(symbolpn, 1);

		index = pppicker.getNextFreeMicrograph();
		if (index == -1)
			index = 0;
		untiltedmic = pppicker.getMicrographs().get(index);

		colorbt.addActionListener(new ColorActionListener());

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
							pppicker.setColor(colorChooser.getColor());
						}
					}, // OK button handler
					null); // no CANCEL button handler
			XmippIJUtil.setLocation(position, 0.5f, dialog);
			dialog.setVisible(true);
		}
	}


	private void initMicrographsPane()
	{
		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(0, 5, 0, 5);
		constraints.anchor = GridBagConstraints.NORTHWEST;
		micrographpn = new JPanel(new GridBagLayout());
		micrographpn.setBorder(BorderFactory.createTitledBorder("Micrographs"));
		JScrollPane sp = new JScrollPane();
		micrographsmd = new MicrographPairsTableModel(this);
		micrographstb = new JTable(micrographsmd);
		micrographstb.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		micrographstb.getColumnModel().getColumn(0).setPreferredWidth(35);
		micrographstb.getColumnModel().getColumn(1).setPreferredWidth(120);
		micrographstb.getColumnModel().getColumn(2).setPreferredWidth(120);
		micrographstb.getColumnModel().getColumn(3).setPreferredWidth(60);
		micrographstb.getColumnModel().getColumn(4).setPreferredWidth(60);
		micrographstb.setPreferredScrollableViewportSize(new Dimension(395, 304));
		micrographstb.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

		sp.setViewportView(micrographstb);
		micrographpn.add(sp, XmippIJUtil.getConstraints(constraints, 0, 0, 1));
		JPanel infopn = new JPanel();
		upslb = new JLabel(Integer.toString(pppicker.getUntiltedNumber()));
		infopn.add(new JLabel("Particles:"));
		infopn.add(upslb);
		micrographpn.add(infopn, XmippIJUtil.getConstraints(constraints, 0, 1, 1));
		JPanel buttonspn = new JPanel(new FlowLayout(FlowLayout.LEFT));
		resetbt = new JButton("Reset");
		buttonspn.add(resetbt);
		micrographpn.add(buttonspn, XmippIJUtil.getConstraints(constraints, 0, 2, 2));
		resetbt.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				pppicker.resetMicrograph(untiltedmic);
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
				if (micrographstb.getSelectedRow() == -1)
					return;// Probably from fireTableDataChanged raised
				index = micrographstb.getSelectedRow();
				// by me.
				untiltedmic.releaseImage();
				untiltedmic = pppicker.getMicrographs().get(index);
				anglesmi.setEnabled(untiltedmic.getAddedCount() >= 4);
				initializeCanvas();
				saveChanges();
				pack();
			}
		});
		micrographstb.getSelectionModel().setSelectionInterval(index, index);

	}

	protected void saveChanges()
	{
		pppicker.saveData();
		setChanged(false);
	}

	public void switchSize(int size)
	{
		super.switchSize(size);
		tiltedcanvas.repaint();
	}

	public void setChanged(boolean changed)
	{
		pppicker.setChanged(changed);
		savemi.setEnabled(changed);
	}

	public void updateMicrographsModel()
	{
		super.updateMicrographsModel();
		micrographsmd.fireTableRowsUpdated(index, index);
		micrographstb.setRowSelectionInterval(index, index);
		upslb.setText(Integer.toString(pppicker.getUntiltedNumber()));
		anglesmi.setEnabled(untiltedmic.getAddedCount() >= 4);
	}

	void initializeCanvas()
	{
		if (canvas == null)
		{
			canvas = new UntiltedMicrographCanvas(this);
			tiltedcanvas = new TiltedMicrographCanvas(this);
			if(!untiltedmic.getParticles().isEmpty())
				canvas.setActive(untiltedmic.getParticles().get(untiltedmic.getParticles().size() - 1));//needs both canvas to be initialized
			if(!pppicker.getFilters().isEmpty())
				IJ.runMacro(pppicker.getFiltersMacro());
		}
		else
		{
			canvas.updateMicrograph();
			tiltedcanvas.updateMicrograph();
		}
	}

	public ParticlePickerCanvas getCanvas()
	{
		return canvas;
	}

	public UntiltedMicrograph getMicrograph()
	{
		return untiltedmic;
	}

	public int getParticleSize()
	{
		return sizesl.getValue();
	}

	public Color getColor()
	{
		return pppicker.getColor();
	}

	

	@Override
	public Family getFamily()
	{
		return pppicker.getFamily();
	}


	@Override
	public List<? extends TrainingParticle> getParticles()
	{
		return untiltedmic.getParticles();
	}


	public TiltedMicrographCanvas getTiltedCanvas()
	{
		return tiltedcanvas;
	}
	
	public boolean drawAngles()
	{
		return anglesmi.isSelected();
	}
	
	

	@Override
	public void changeShapes()
	{
		canvas.repaint();
		tiltedcanvas.repaint();
		
	}
	
	public  void importParticlesXmipp24(String projectdir)
	{
		super.importParticlesXmipp24(projectdir);
		tiltedcanvas.repaint();
		canvas.setActive(null);
		
	}
	
	protected void importParticlesXmipp30(String path)
	{
		super.importParticlesXmipp30(path);
		untiltedmic.initAligner();
		tiltedcanvas.repaint();
	}
	
	public void importParticlesFromXmipp24Files(String ufile, String tfile)
	{
		getParticlePicker().importParticlesFromXmipp24Files(untiltedmic, ufile, tfile, true);
		untiltedmic.initAligner();
		setChanged(true);
		getCanvas().repaint();
		updateMicrographsModel();
		tiltedcanvas.repaint();
		canvas.setActive(null);
		
	}



	@Override
	protected void displayImportDialog()
	{
		new ImportParticlesFromFilesJDialog(this, true);
		
	}
	
	


}
