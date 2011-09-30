package trainingpicker.gui;

import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.Panel;
import java.awt.ScrollPane;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.util.ArrayList;
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

public class ParticlesJDialog extends JDialog
{

	private ParticlePickerJFrame parent;
	private List<TrainingParticle> particles;
	private int rows, columns;
	private int side;
	private List<ParticleCanvas> canvass;

	public ParticlesJDialog(ParticlePickerJFrame parent, TrainingMicrograph micrograph)
	{
		super(parent);
		this.parent = parent;
		particles = micrograph.getFamilyData(parent.getFamily()).getParticles();
		side = (int)(parent.getFamily().getSize() * parent.getMagnification());
		columns = 800 / side;
		rows = particles.size() / columns + 1;
		canvass = new ArrayList<ParticleCanvas>();
		initComponents();
	}

	private void initComponents()
	{
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Particles");
		setResizable(false);
		GridBagConstraints constraints = new GridBagConstraints();
		ScrollPane sp = new ScrollPane();
		

		Panel particlespn = new Panel(new GridBagLayout());
		int index = 0;

		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns; j++, index++)
			{
				if (index == particles.size())
					break;
				canvass.add(new ParticleCanvas(particles.get(index), parent));
				particlespn.add(canvass.get(index), WindowUtils.getConstraints(constraints, j, i, 1));
			}
		int width = side * columns;
		int height = side * 2;
		sp.setPreferredSize(new Dimension(width + 20, height));
		sp.add(particlespn);
		add(sp);
		pack();
		setVisible(true);

	}

	

}
