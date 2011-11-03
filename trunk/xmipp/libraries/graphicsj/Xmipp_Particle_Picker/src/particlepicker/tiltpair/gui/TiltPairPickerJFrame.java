package particlepicker.tiltpair.gui;

import ij.IJ;
import ij.ImageJ;

import java.awt.Color;
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
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JColorChooser;
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
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import particlepicker.Family;
import particlepicker.ParticlePicker;
import particlepicker.ParticlePickerCanvas;
import particlepicker.ParticlePickerJFrame;
import particlepicker.Shape;
import particlepicker.Tool;
import particlepicker.WindowUtils;
import particlepicker.tiltpair.model.TiltPairPicker;
import particlepicker.tiltpair.model.UntiltedMicrograph;
import particlepicker.training.gui.ColorIcon;
import particlepicker.training.model.TrainingParticle;
import particlepicker.training.model.TrainingPicker;




public class TiltPairPickerJFrame extends ParticlePickerJFrame 
{

	private TiltPairPicker pppicker;
	private JSlider sizesl;
	private JFormattedTextField sizetf;
	private JMenuBar mb;
	private Color color;
	private JPanel particlespn;
	private JPanel micrographpn;
	private JTable micrographstb;
	private UntiltedMicrographCanvas canvas;
	private MicrographPairsTableModel micrographsmd;
	private UntiltedMicrograph untiltedmic;
	private JButton colorbt;
	private double position;
	private int index;
	private JButton resetbt;
	private JLabel upslb;
	private TiltedMicrographCanvas tiltedcanvas;
	private JCheckBoxMenuItem anglesmi;

	public TiltPairPicker getParticlePicker()
	{
		return pppicker;
	}
	


	public TiltPairPickerJFrame(TiltPairPicker pppicker)
	{

		this.pppicker = pppicker;
		initComponents();
		// initializeCanvas();

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
		add(particlespn, WindowUtils.getConstraints(constraints, 0, 1, 3));
		initMicrographsPane();
		add(micrographpn, WindowUtils.getConstraints(constraints, 0, 2, 3));

		pack();
		position = 0.95;
		WindowUtils.centerScreen(position, 0.5, this);
		setVisible(true);
	}
	


	public void initMenuBar()
	{
		mb = new JMenuBar();

		// Setting menus
		JMenu filemn = new JMenu("File");
		
		JMenu viewmn = new JMenu("View");
		JMenu helpmn = new JMenu("Help");
		mb.add(filemn);
		mb.add(filtersmn);
		mb.add(viewmn);
		mb.add(helpmn);

		// Setting menu items
		
		savemi.setEnabled(pppicker.isChanged());
		filemn.add(savemi);
		anglesmi = new JCheckBoxMenuItem("Angles");
		anglesmi.setEnabled(false);
		anglesmi.addActionListener(new ActionListener()
		{
			
			@Override
			public void actionPerformed(ActionEvent e)
			{
				
				
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
		int size = pppicker.getSize();
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

		particlespn.add(fieldspn, 0);
		initSymbolPane();
		particlespn.add(symbolpn, 1);

		index = pppicker.getNextFreeMicrograph();
		if (index == -1)
			index = 0;
		untiltedmic = pppicker.getMicrographs().get(index);

		colorbt.addActionListener(new ColorActionListener());

		sizetf.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				int size = ((Number) sizetf.getValue()).intValue();
				switchSize(size);
			}
		});

		sizesl.addChangeListener(new ChangeListener()
		{

			@Override
			public void stateChanged(ChangeEvent e)
			{
				int size = sizesl.getValue();
				switchSize(size);
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
							pppicker.setColor(colorChooser.getColor());
						}
					}, // OK button handler
					null); // no CANCEL button handler
			WindowUtils.centerScreen(position, 0.5, dialog);
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
		micrographstb.setPreferredScrollableViewportSize(new Dimension(335, 304));
		micrographstb.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

		sp.setViewportView(micrographstb);
		micrographpn.add(sp, WindowUtils.getConstraints(constraints, 0, 0, 1));
		JPanel infopn = new JPanel();
		upslb = new JLabel(Integer.toString(pppicker.getUntiltedNumber()));
		infopn.add(new JLabel("Particles:"));
		infopn.add(upslb);
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

	void switchSize(int size)
	{
		sizetf.setText(Integer.toString(size));
		sizesl.setValue(size);
		canvas.repaint();
		tiltedcanvas.repaint();
		pppicker.setSize(size);
		setChanged(true);
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
	}

	void initializeCanvas()
	{
		if (canvas == null)
		{
			canvas = new UntiltedMicrographCanvas(this);
			tiltedcanvas = new TiltedMicrographCanvas(this);
			if(!untiltedmic.getParticles().isEmpty())
				canvas.setActive(untiltedmic.getParticles().get(untiltedmic.getParticles().size() - 1));//needs both canvas to be initialized
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
}
