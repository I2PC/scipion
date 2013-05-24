package xmipp.viewer.particlepicker.tiltpair.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;

import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.Family;
import xmipp.viewer.particlepicker.Format;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.ParticlesJDialog;
import xmipp.viewer.particlepicker.tiltpair.model.TiltPairPicker;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedParticle;
import xmipp.viewer.particlepicker.training.model.FamilyState;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;




public class TiltPairPickerJFrame extends ParticlePickerJFrame 
{

	private TiltPairPicker pppicker;
	private JMenuBar mb;
	private JPanel particlespn;
	private JPanel micrographpn;
	private UntiltedMicrographCanvas canvas;
	private MicrographPairsTableModel micrographsmd;
	
	
	private float position;
	private int index;
	private JLabel upslb;
	private TiltedMicrographCanvas tiltedcanvas;
	private JCheckBoxMenuItem anglesmi;
	private JMenuItem importffilesmi;
	private ImportParticlesFromFilesTiltPairJDialog importffilesjd;

	public TiltPairPicker getParticlePicker()
	{
		return pppicker;
	}

	public TiltPairPickerJFrame(TiltPairPicker picker)
	{
		super(picker);
		pppicker = picker;
		initComponents();
		enableEdition(picker.getMode() != FamilyState.ReadOnly);
		setChanged(false);
	}

	private void initComponents()
	{
		setResizable(false);
		setTitle("Xmipp Particle Pair Picker");
		initMenuBar();
		setJMenuBar(mb);

		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(0, 5, 0, 5);
		constraints.anchor = GridBagConstraints.WEST;
		setLayout(new GridBagLayout());

		initParticlesPane();
		add(particlespn, XmippWindowUtil.getConstraints(constraints, 0, 1, 3));
		initMicrographsPane();
		add(micrographpn, XmippWindowUtil.getConstraints(constraints, 0, 2, 3));
		JPanel actionspn = new JPanel(new FlowLayout(FlowLayout.RIGHT));
		actionspn.add(savebt);
		actionspn.add(saveandexitbt);
		add(actionspn, XmippWindowUtil.getConstraints(constraints, 0, 3, 3, GridBagConstraints.HORIZONTAL));

		pack();
		position = 0.95f;
		XmippWindowUtil.setLocation(position, 0.5f, this);
		setVisible(true);
	}
	


	public void initMenuBar()
	{
		mb = new JMenuBar();

		filemn.add(importffmi);
		if (pppicker.getFamily().getStep() != FamilyState.Manual)
			importffmi.setEnabled(false);
		importffilesmi = new JMenuItem("Import Particles From Micrograph");
		importffilesmi.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent arg0) {
				showImportFromFilesDialog();
				
			}
		});
		filemn.add(importffilesmi, 1);
		// Setting menus
		JMenu viewmn = new JMenu("View");
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
	}
	
	protected void showImportFromFilesDialog(){
		if (importffilesjd == null)
			importffilesjd = new ImportParticlesFromFilesTiltPairJDialog(this);
		importffilesjd.showDialog();
	}

	private void initParticlesPane()
	{
		particlespn = new JPanel();
		GridLayout gl = new GridLayout(2, 1);
		particlespn.setLayout(gl);

		particlespn.setBorder(BorderFactory.createTitledBorder("Particles"));

		JPanel fieldspn = new JPanel(new FlowLayout(FlowLayout.LEFT));

		// Setting color
		initColorPane(pppicker.getFamily().getColor());
		fieldspn.add(colorbt);

		// Setting slider
		initSizePane();
		fieldspn.add(sizepn);

		particlespn.add(fieldspn, 0);
		initImagePane();
		particlespn.add(imagepn, 1);

		index = pppicker.getMicrographIndex();

		colorbt.addActionListener(new ColorActionListener());

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
		micrographstb.setModel(micrographsmd);
		formatMicrographsTable();

		sp.setViewportView(micrographstb);
		micrographpn.add(sp, XmippWindowUtil.getConstraints(constraints, 0, 0, 1));
		JPanel infopn = new JPanel();
		upslb = new JLabel(Integer.toString(pppicker.getUntiltedNumber()));
		infopn.add(new JLabel("Particles:"));
		infopn.add(upslb);
		micrographpn.add(infopn, XmippWindowUtil.getConstraints(constraints, 0, 1, 1));
		JPanel buttonspn = new JPanel(new FlowLayout(FlowLayout.LEFT));
		
		buttonspn.add(resetbt);
		micrographpn.add(buttonspn, XmippWindowUtil.getConstraints(constraints, 0, 2, 2));
		

	}
	
	private void formatMicrographsTable() {
		micrographstb.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		micrographstb.getColumnModel().getColumn(0).setPreferredWidth(35);
		micrographstb.getColumnModel().getColumn(1).setPreferredWidth(220);
		micrographstb.getColumnModel().getColumn(2).setPreferredWidth(220);
		micrographstb.getColumnModel().getColumn(3).setPreferredWidth(60);
		micrographstb.getColumnModel().getColumn(4).setPreferredWidth(60);
		micrographstb.setPreferredScrollableViewportSize(new Dimension(595, 304));
		micrographstb.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		if(index != -1)
			micrographstb.setRowSelectionInterval(index, index);
	}

	
	

	public void updateSize(int size)
	{
		super.updateSize(size);
		tiltedcanvas.repaint();
	}

	public void setChanged(boolean changed)
	{
		pppicker.setChanged(changed);
		savemi.setEnabled(changed);
		savebt.setEnabled(changed);
	}

	public void updateMicrographsModel(boolean all)
	{
		if (particlesdialog != null)
			loadParticles();
		if(all)
			micrographsmd.fireTableRowsUpdated(0, micrographsmd.getRowCount() - 1 );
		else
			micrographsmd.fireTableRowsUpdated(index, index);
		micrographstb.setRowSelectionInterval(index, index);
		upslb.setText(Integer.toString(pppicker.getUntiltedNumber()));
		anglesmi.setEnabled(pppicker.getMicrograph().getAddedCount() >= 4);
	}

	void initializeCanvas()
	{
		if (canvas == null)
		{
			canvas = new UntiltedMicrographCanvas(this);
			tiltedcanvas = new TiltedMicrographCanvas(this);
			List<UntiltedParticle> particles = getMicrograph().getParticles();
			if(!particles.isEmpty())
				canvas.refreshActive(particles.get(particles.size() - 1));//needs both canvas to be initialized
		}
		else
		{
			canvas.updateMicrograph();
			tiltedcanvas.updateMicrograph();
						
		}
		canvas.display();
		tiltedcanvas.display(0.7f, 0);
		updateZoom();
		
		if(usezoombt.isSelected())
			tiltedcanvas.setZoom(getZoom());
	}

	public ParticlePickerCanvas getCanvas()
	{
		return canvas;
	}

	public UntiltedMicrograph getMicrograph()
	{
		return pppicker.getMicrograph();
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
	public List<? extends TrainingParticle> getAvailableParticles()
	{
		return getMicrograph().getParticles();
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

	public String importParticlesFromFolder(Format format, String dir, float scale, boolean invertx, boolean inverty)
	{
		String result = pppicker.importParticlesFromFolder(dir, format, scale, invertx, inverty);
		getCanvas().repaint();
		updateMicrographsModel(true);
		getCanvas().refreshActive(null);
		tiltedcanvas.repaint();
		return result;
	}
	


	@Override
	protected void resetMicrograph()
	{
		pppicker.resetMicrograph(getMicrograph());
		canvas.refreshActive(null);
		updateMicrographsModel();
		setChanged(true);
		
	}
	
	@Override
	protected void reloadImage()
	{
		getCanvas().getMicrograph().releaseImage();
		getCanvas().updateMicrograph();
		getTiltedCanvas().getMicrograph().releaseImage();
		getTiltedCanvas().updateMicrograph();
		canvas.display();
		getTiltedCanvas().display();
		
	}

	@Override
	protected void loadMicrograph()
	{
		if (this.micrographstb.getSelectedRow() == -1)
			return;// Probably from fireTableDataChanged raised
		if(index == micrographstb.getSelectedRow() && canvas != null && canvas.getIw().isVisible())//micrograph open, no need to reopen
			return;
		pppicker.saveData(getMicrograph());
		getMicrograph().releaseImage();
		
		index = micrographstb.getSelectedRow();
		pppicker.setMicrograph(pppicker.getMicrographs().get(index));
		pppicker.saveConfig();
		setChanged(false);
		anglesmi.setEnabled(getMicrograph().getAddedCount() >= 4);
		initializeCanvas();
		
		pack();
		if (particlesdialog != null)
			loadParticles();
	}

	
	@Override
	protected void openHelpURl() {
		XmippWindowUtil.openURI("http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Micrograph__picking_v3");
	}


	
	public void importParticlesFromFiles(Format format, String file1, String file2, float scale, boolean invertx, boolean inverty){
			
			getMicrograph().reset();
			String result = pppicker.importParticlesFromFiles(file1, file2, format, getMicrograph(), scale, invertx, inverty);
			pppicker.saveData(getMicrograph());
			setChanged(false);
			getCanvas().repaint();
			tiltedcanvas.repaint();
			updateMicrographsModel();
			updateSize(getFamily().getSize());
			canvas.refreshActive(null);
	}
	
	@Override
	public String importParticles(Format format, String dir, float scale, boolean invertx, boolean inverty)
	{
		return importParticlesFromFolder(format, dir, scale, invertx, inverty);
		
	}

	@Override
	public ParticlesJDialog initParticlesJDialog()
	{
		return new TiltPairParticlesJDialog(this);
	}

}//class TiltPairPickerJFrame
