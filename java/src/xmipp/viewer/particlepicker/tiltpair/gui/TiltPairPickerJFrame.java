package xmipp.viewer.particlepicker.tiltpair.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
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

import xmipp.utils.XmippDialog;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.Format;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.ParticlesJDialog;
import xmipp.viewer.particlepicker.tiltpair.model.TiltPairPicker;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedParticle;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.ManualParticle;

public class TiltPairPickerJFrame extends ParticlePickerJFrame {

	private TiltPairPicker tppicker;
	private JMenuBar mb;
	private JPanel micrographpn;
	private UntiltedMicrographCanvas canvas;
	private MicrographPairsTableModel micrographsmd;

	private float position;
	
	private JLabel upslb;
	private TiltedMicrographCanvas tiltedcanvas;
	private JCheckBoxMenuItem anglesmi;
	private JMenuItem importffilesmi;
	private ImportParticlesFromFilesTiltPairJDialog importffilesjd;

	public TiltPairPicker getParticlePicker() {
		return tppicker;
	}

	public TiltPairPickerJFrame(TiltPairPicker picker) {
		super(picker);
		tppicker = picker;
		initComponents();

		setChanged(false);
	}

	private void initComponents() {
		setResizable(false);
		setTitle("Xmipp Tilt Pair Picker - " + tppicker.getMode());
		initMenuBar();
		setJMenuBar(mb);

		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(0, 5, 0, 5);
		constraints.anchor = GridBagConstraints.WEST;
		setLayout(new GridBagLayout());

		initToolBar();
		add(tb);
		initShapePane();
		JPanel shapepn2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
		shapepn2.add(new JLabel("Shape:"));
		shapepn2.add(shapepn);
		add(shapepn2, XmippWindowUtil.getConstraints(constraints, 0, 1, 3));
		initMicrographsPane();
		add(micrographpn, XmippWindowUtil.getConstraints(constraints, 0, 2, 3));
		JPanel actionspn = new JPanel(new FlowLayout(FlowLayout.RIGHT));
		actionspn.add(savebt);
		actionspn.add(saveandexitbt);
		add(actionspn, XmippWindowUtil.getConstraints(constraints, 0, 3, 3,
				GridBagConstraints.HORIZONTAL));

		actionspn.add(savebt);
		actionspn.add(saveandexitbt);
		add(actionspn, XmippWindowUtil.getConstraints(constraints, 0, 3, 3, 1,
				GridBagConstraints.HORIZONTAL));

		enableEdition(tppicker.getMode() != Mode.ReadOnly);
		pack();
		position = 0.95f;
		XmippWindowUtil.setLocation(position, 0.5f, this);
		setVisible(true);
	}

	protected void enableEdition(boolean enable) {
		super.enableEdition(enable);
		importffilesmi.setEnabled(enable);
	}

	public void initMenuBar() {
		mb = new JMenuBar();

		filemn.add(importmi);
		if (tppicker.getMode() != Mode.Manual)
			importmi.setEnabled(false);
		importffilesmi = new JMenuItem("Import Micrograph Particles");
		importffilesmi.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent arg0) {
				showImportFromFilesDialog();

			}
		});
		filemn.add(importffilesmi, 1);
		filemn.add(exitmi);
		// Setting menus
		JMenu viewmn = new JMenu("View");
		mb.add(filemn);
		mb.add(filtersmn);
		mb.add(viewmn);
		mb.add(helpmn);

		anglesmi = new JCheckBoxMenuItem("Angles");
		// anglesmi.setEnabled(false);
		anglesmi.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				if (anglesmi.isSelected()
						&& !tppicker.getMicrograph().anglesAvailable())
					XmippDialog
							.showInfo(TiltPairPickerJFrame.this,
									"Angles will not be accurate until more than 10 pairs are registered");
				canvas.repaint();
				tiltedcanvas.repaint();
			}
		});
		viewmn.add(anglesmi);
		viewmn.add(pmi);
		viewmn.add(ijmi);
	}

	protected void showImportFromFilesDialog() {
		if (importffilesjd == null)
			importffilesjd = new ImportParticlesFromFilesTiltPairJDialog(this);
		importffilesjd.showDialog();
	}

	private void initMicrographsPane() {
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
		micrographpn.add(sp,
				XmippWindowUtil.getConstraints(constraints, 0, 0, 1));
		JPanel infopn = new JPanel();
		upslb = new JLabel(Integer.toString(tppicker.getUntiltedNumber()));
		infopn.add(new JLabel("Particles:"));
		infopn.add(upslb);
		micrographpn.add(infopn,
				XmippWindowUtil.getConstraints(constraints, 0, 1, 1));
		JPanel buttonspn = new JPanel(new FlowLayout(FlowLayout.LEFT));

		buttonspn.add(resetbt);
		micrographpn.add(buttonspn,
				XmippWindowUtil.getConstraints(constraints, 0, 2, 2));

	}

	private void formatMicrographsTable() {
		micrographstb.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		micrographstb.getColumnModel().getColumn(0).setPreferredWidth(35);
		micrographstb.getColumnModel().getColumn(1).setPreferredWidth(220);
		micrographstb.getColumnModel().getColumn(2).setPreferredWidth(220);
		micrographstb.getColumnModel().getColumn(3).setPreferredWidth(60);
		micrographstb.getColumnModel().getColumn(4).setPreferredWidth(60);
		micrographstb
				.setPreferredScrollableViewportSize(new Dimension(595, 304));
		micrographstb.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		int index = tppicker.getMicrographIndex();
		if (index != -1)
			micrographstb.setRowSelectionInterval(index, index);
	}

	public void updateSize(int size) {
		super.updateSize(size);
		tiltedcanvas.repaint();
	}

	public void setChanged(boolean changed) {
		tppicker.setChanged(changed);
		savemi.setEnabled(changed);
		savebt.setEnabled(changed);
	}

	public void updateMicrographsModel(boolean all) {
		if (particlesdialog != null)
			loadParticles();
		int index = tppicker.getMicrographIndex();
		if (all)
			micrographsmd.fireTableRowsUpdated(0,
					micrographsmd.getRowCount() - 1);
		else
			micrographsmd.fireTableRowsUpdated(index, index);
		micrographstb.setRowSelectionInterval(index, index);
		upslb.setText(Integer.toString(tppicker.getUntiltedNumber()));
		// anglesmi.setEnabled(tppicker.getMicrograph().getAddedCount() >= 4);
	}

	void initializeCanvas() {
		if (canvas == null) {
			canvas = new UntiltedMicrographCanvas(this);
			tiltedcanvas = new TiltedMicrographCanvas(this);
			List<UntiltedParticle> particles = getMicrograph().getParticles();
			if (!particles.isEmpty())
				canvas.refreshActive(particles.get(particles.size() - 1));// needs
																			// both
																			// canvas
																			// to
																			// be
																			// initialized
		} else {
			canvas.updateMicrograph();
			tiltedcanvas.updateMicrograph();

		}
		canvas.display();

		tiltedcanvas.display();
		updateZoom();
		tiltedcanvas.getIw().setLocation(canvas.getIw().getWidth(), 0);
		if (usezoombt.isSelected())
			tiltedcanvas.setZoom(getZoom());
	}

	public ParticlePickerCanvas getCanvas() {
		return canvas;
	}

	public UntiltedMicrograph getMicrograph() {
		return tppicker.getMicrograph();
	}

	public int getParticleSize() {
		return sizesl.getValue();
	}

	public Color getColor() {
		return tppicker.getColor();
	}

	@Override
	public List<? extends ManualParticle> getAvailableParticles() {
		return getMicrograph().getParticles();
	}

	public TiltedMicrographCanvas getTiltedCanvas() {
		return tiltedcanvas;
	}

	public boolean drawAngles() {
		return anglesmi.isSelected();
	}

	@Override
	public void changeShapes() {
		canvas.repaint();
		tiltedcanvas.repaint();
	}

	public String importParticlesFromFolder(Format format, String dir,
			float scale, boolean invertx, boolean inverty) {
		String result = tppicker.importParticlesFromFolder(dir, format, scale,
				invertx, inverty);
		getCanvas().repaint();
		updateMicrographsModel(true);
		getCanvas().refreshActive(null);
		tiltedcanvas.repaint();
		return result;
	}

	@Override
	protected void resetMicrograph() {
		tppicker.resetMicrograph(getMicrograph());
		canvas.refreshActive(null);
		updateMicrographsModel();

	}

	@Override
	protected void reloadImage() {
		getCanvas().getMicrograph().releaseImage();
		getCanvas().updateMicrograph();
		getTiltedCanvas().getMicrograph().releaseImage();
		getTiltedCanvas().updateMicrograph();
		canvas.display();
		getTiltedCanvas().display();

	}

	@Override
	protected void loadMicrograph() {
		if (this.micrographstb.getSelectedRow() == -1)
			return;// Probably from fireTableDataChanged raised
		int index = tppicker.getMicrographIndex();
		if (index == micrographstb.getSelectedRow() && canvas != null
				&& canvas.getIw().isVisible())// micrograph open, no need to
												// reopen
			return;
		if (tppicker.isChanged()) 
			tppicker.saveData(getMicrograph());
		getMicrograph().releaseImage();

		index = micrographstb.getSelectedRow();
		tppicker.setMicrograph(tppicker.getMicrographs().get(index));
		tppicker.saveConfig();
		
		
		// anglesmi.setEnabled(getMicrograph().getAddedCount() >= 4);
		initializeCanvas();

		pack();
		if (particlesdialog != null)
			loadParticles();
	}

	@Override
	protected void openHelpURl() {
		XmippWindowUtil
				.openURI("http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Micrograph_tiltpair_picking_v3");
	}

	public String importParticlesFromFiles(Format format, String file1,
			String file2, float scale, boolean invertx, boolean inverty) {

		getMicrograph().reset(tppicker);
		String result = tppicker.importParticlesFromFiles(file1, file2, format,
				getMicrograph(), scale, invertx, inverty);
		tppicker.saveData(getMicrograph());
		getCanvas().repaint();
		tiltedcanvas.repaint();
		updateMicrographsModel();
		updateSize(tppicker.getSize());
		canvas.refreshActive(null);
		return result;
	}

	@Override
	public String importParticles(Format format, String dir, float scale,
			boolean invertx, boolean inverty) {
		return importParticlesFromFolder(format, dir, scale, invertx, inverty);

	}

	@Override
	public ParticlesJDialog initParticlesJDialog() {
		return new TiltPairParticlesJDialog(this);
	}

}// class TiltPairPickerJFrame
