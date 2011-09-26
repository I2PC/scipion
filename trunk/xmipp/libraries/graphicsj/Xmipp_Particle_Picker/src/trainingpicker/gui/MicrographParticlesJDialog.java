package trainingpicker.gui;

import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.util.List;

import javax.swing.AbstractListModel;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;

import trainingpicker.model.TrainingMicrograph;
import trainingpicker.model.TrainingParticle;
import xmipp.Particle;

public class MicrographParticlesJDialog extends JDialog
{

	private ParticlePickerJFrame parent;
	private TrainingMicrograph micrograph;
	private int size;
	private List<TrainingParticle> particles;
	private int width;

	public MicrographParticlesJDialog(ParticlePickerJFrame parent, TrainingMicrograph micrograph)
	{
		super(parent);
		this.parent = parent;
		this.micrograph = micrograph;
		particles = micrograph.getFamilyData(parent.getFamily()).getManualParticles();
		this.size = parent.getFamily().getSize();
		width = (int) Math.round(Math.sqrt(particles.size()) + 1);
		initComponents();

	}

	private void initComponents()
	{
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Cells Info");
		GridBagConstraints constraints = new GridBagConstraints();
		setLayout(new GridBagLayout());
		ImagePlus imp;
		int index = 0;

		for (int j = 0; j < width; j++)
			for (int i = 0; i < width; i++, index++)
			{
				if (index == particles.size())
					break;
				imp = particles.get(index).getImage();
				add(new ParticleCanvas(imp), WindowUtils.getConstraints(constraints, i, j, 1));
			}
		pack();
		centerScreen();
		setVisible(true);
	}

	private void centerScreen()
	{
		// TODO Auto-generated method stub

	}

}
