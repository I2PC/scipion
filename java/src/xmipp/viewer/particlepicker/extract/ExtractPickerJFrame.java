package xmipp.viewer.particlepicker.extract;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.border.TitledBorder;
import javax.swing.table.AbstractTableModel;

import xmipp.ij.commons.Tool;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.ColorHelper;
import xmipp.viewer.particlepicker.Format;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.ParticlesJDialog;
import xmipp.viewer.particlepicker.PickerParticle;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.windows.GalleryJFrame;
import xmipp.viewer.windows.ImagesWindowFactory;

public class ExtractPickerJFrame extends ParticlePickerJFrame
{

	private ExtractParticlePicker picker;
	private JMenuBar mb;
	private JPanel micrographpn;
	private JButton iconbt;
	private ExtractMicrographsTableModel micrographsmd;
	private JLabel particleslb;
	private int index;
	private ExtractCanvas canvas;
	private ColorHelper[] colorby;
	private JComboBox scorescb;
	private GalleryJFrame galleryfr;
	private JLabel minlb;
	private JLabel maxlb;
	private ExtractParticle active;
	private JPanel colorlegendpn;

	public ExtractPickerJFrame(ParticlePicker picker, GalleryJFrame galleryfr)
	{
		super(picker);
		this.picker = (ExtractParticlePicker)picker;
		this.galleryfr = galleryfr;
		initComponents();
		// TODO Auto-generated constructor stub
	}

	private void initComponents()
	{
		setResizable(false);
		setTitle("Xmipp Particle Picker - " + picker.getMode());
		initMenuBar();
		setJMenuBar(mb);
		savemi.setEnabled(false);

		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(0, 5, 0, 5);
		constraints.anchor = GridBagConstraints.WEST;
		setLayout(new GridBagLayout());

		initToolBar();
		add(tb, XmippWindowUtil.getConstraints(constraints, 0, 1, 2, 1, GridBagConstraints.NONE));
		add(new JLabel("Shape:"), XmippWindowUtil.getConstraints(constraints, 0, 2));
		initShapePane();
		add(shapepn, XmippWindowUtil.getConstraints(constraints, 1, 2));
		add(new JLabel("Color By:"), XmippWindowUtil.getConstraints(constraints, 0, 3));
		initColorLegendPane();
		add(colorlegendpn, XmippWindowUtil.getConstraints(constraints, 1, 3));
		initMicrographsPane();

		initMicrographsPane();
		add(micrographpn, XmippWindowUtil.getConstraints(constraints, 0, 4, 3));

		pack();
		float positionx = 0.995f;
		XmippWindowUtil.setLocation(positionx, 0.25f, this);
		setVisible(true);
	}
	
	
	public boolean isPickingAvailable(MouseEvent e)
	{
		return false;
	}
	
	private void initColorLegendPane()
	{

		colorlegendpn = new JPanel(new FlowLayout(FlowLayout.LEFT));

		colorlegendpn.add(new JLabel("Color by:"));
		colorby = picker.getColumns();
		scorescb = new JComboBox(colorby);
		scorescb.addActionListener(new ActionListener()
		{
			
			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				getCanvas().repaint();
				minlb.setText(String.format("%.2f", getColorHelper().getMin()));
				maxlb.setText(String.format("%.2f", getColorHelper().getMax()));
			}
		});
		colorlegendpn.add(scorescb);
		minlb = new JLabel(String.format("%.2f", getColorHelper().getMin()));
		maxlb = new JLabel(String.format("%.2f", getColorHelper().getMax()));
		
		
		colorlegendpn.add(minlb);
		colorlegendpn.add(ColorHelper.getColorMap());
		colorlegendpn.add(maxlb);
		
		
		
		index = picker.getMicrographIndex();

	}
	
	ColorHelper getColorHelper()
	{
		return (ColorHelper)scorescb.getSelectedItem();
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
					ImagesWindowFactory.openCTFWindow(getMicrograph().getPSDImage(), getMicrograph().getCTF(), getMicrograph().getPSD());

			}
		});
		ctfpn.add(iconbt);
		micrographsmd = new ExtractMicrographsTableModel();
		micrographstb.setModel(micrographsmd);
		micrographstb.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		micrographstb.getColumnModel().getColumn(0).setPreferredWidth(35);
		micrographstb.getColumnModel().getColumn(1).setPreferredWidth(245);
		micrographstb.getColumnModel().getColumn(2).setPreferredWidth(70);
		micrographstb.getColumnModel().getColumn(3).setPreferredWidth(70);
		micrographstb.setPreferredScrollableViewportSize(new Dimension(420, 304));
		micrographstb.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		if (index != -1)
			micrographstb.setRowSelectionInterval(index, index);

		sp.setViewportView(micrographstb);
		micrographpn.add(sp, XmippWindowUtil.getConstraints(constraints, 0, 0, 1));
		micrographpn.add(ctfpn, XmippWindowUtil.getConstraints(constraints, 1, 0, 1));
		JPanel infopn = new JPanel();
		particleslb = new JLabel(Integer.toString(picker.getParticlesTotal()));
		
		infopn.add(new JLabel("Particles:"));
		infopn.add(particleslb);
		micrographpn.add(infopn, XmippWindowUtil.getConstraints(constraints, 0, 1, 1));
		JPanel buttonspn = new JPanel(new FlowLayout(FlowLayout.LEFT));

		buttonspn.add(resetbt);
		micrographpn.add(buttonspn, XmippWindowUtil.getConstraints(constraints, 0, 2, 2));

		
	}

	private void initMenuBar()
	{
		mb = new JMenuBar();

		
		JMenu windowmn = new JMenu("Window");
		JMenu helpmn = new JMenu("Help");
		mb.add(filemn);
		mb.add(filtersmn);
		mb.add(windowmn);
		mb.add(helpmn);
		// importffilemi.setText("Import from File...");

		windowmn.add(pmi);
		windowmn.add(ijmi);
		helpmn.add(hcontentsmi);
	
	}


	@Override
	protected void loadMicrograph()
	{
		if (micrographstb.getSelectedRow() == -1)
			return;// Probably from fireTableDataChanged raised
		// is same micrograph??
		if (index == micrographstb.getSelectedRow() && canvas != null && canvas.getIw().isVisible())
			return;
		if (picker.isChanged())
			picker.saveData(getMicrograph());// Saving changes when switching
		index = micrographstb.getSelectedRow();
		picker.getMicrograph().releaseImage();
		picker.setMicrograph(picker.getMicrographs().get(index));
		setChanged(false);
		initializeCanvas();
		iconbt.setIcon(picker.getMicrograph().getCTFIcon());
		
		pack();
	}

	private void initializeCanvas()
	{
		if (canvas == null)
			canvas = new ExtractCanvas(this);
		else
			canvas.updateMicrograph();
		
		canvas.display();
		updateZoom();
		
	}

	@Override
	protected void openHelpURl()
	{
		// TODO Auto-generated method stub

	}

	@Override
	protected void resetMicrograph()
	{
		picker.getMicrograph().getParticles().clear();
		canvas.refreshActive(null);
		updateMicrographsModel();

	}

	@Override
	protected void reloadImage()
	{
		// TODO Auto-generated method stub

	}

	

	@Override
	public ParticlePickerCanvas getCanvas()
	{
		return canvas;
	}

	@Override
	public void updateMicrographsModel(boolean all)
	{
		if (particlesdialog != null)
			loadParticles();

		if (all)
			micrographsmd.fireTableRowsUpdated(0, micrographsmd.getRowCount() - 1);
		else
			micrographsmd.fireTableRowsUpdated(index, index);

		micrographstb.setRowSelectionInterval(index, index);
		particleslb.setText(Integer.toString(picker.getParticlesTotal()));

	}

	@Override
	public ExtractMicrograph getMicrograph()
	{
		return picker.getMicrograph();
	}

	@Override
	public List<? extends PickerParticle> getAvailableParticles()
	{
		return getMicrograph().getParticles();
	}

	

	@Override
	public ExtractParticlePicker getParticlePicker()
	{
		return picker;
	}

	@Override
	public void setChanged(boolean changed)
	{
		picker.setChanged(changed);
	}

	

	@Override
	public String importParticles(Format format, String dir, float scale, boolean invertx, boolean inverty)
	{
		throw new UnsupportedOperationException();

	}
	
	class ExtractMicrographsTableModel extends AbstractTableModel
	{
		private String[] columns = new String[]{"", "Name", "Particles", "State"};

		@Override
		public int getColumnCount()
		{
			return columns.length;
		}
		
		@Override
		public String getColumnName(int c)
		{
			return columns[c];
		}


		@Override
		public int getRowCount()
		{
			return picker.getMicrographs().size();
		}

		@Override
		public Object getValueAt(int row, int column)
		{
			ExtractMicrograph m = (ExtractMicrograph)picker.getMicrographs().get(row);
			if(column == 0)
				return row + 1;
			if(column == 1)
				return m.getName();
			if(column == 2)
				return m.getParticles().size();
			if(column == 3)
				return Mode.Extract;

			return null;
			
		}
		
	}
	@Override
	public ParticlesJDialog initParticlesJDialog()
	{
		return new ParticlesJDialog(this);
	}

	public void refreshActive(long id, boolean b)
	{
		int index = 0;
		active = null;
		List<ExtractMicrograph> micrographs = picker.getMicrographs();
		for(int i = 0; i < micrographs.size(); i ++)
			for(ExtractParticle p: micrographs.get(i).getParticles())
				if(p.getId() == id)
				{
					index = i;
					active = p;
					active.setEnabled(b);
					break;
				
				}
		if(index != this.index)
			micrographstb.setRowSelectionInterval(index, index);
		if(active != null)
		{
			canvas.refreshActive(active);
			canvas.moveTo(active);
			active = null;
		}

	}

	public void refreshActiveOnGallery(ExtractParticle active)
	{
		if(galleryfr != null)
			galleryfr.refreshActive(active.getId(), active.isEnabled());
		
	}

}
