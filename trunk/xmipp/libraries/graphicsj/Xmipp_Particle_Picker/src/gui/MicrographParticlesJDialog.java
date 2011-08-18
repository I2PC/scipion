package gui;

import ij.ImagePlus;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.LayoutManager;
import java.util.List;

import javax.swing.AbstractListModel;
import javax.swing.JDialog;
import javax.swing.JList;
import javax.swing.JScrollPane;

import model.Micrograph;
import model.Particle;



public class MicrographParticlesJDialog extends JDialog {
	
	private XmippParticlePickerJFrame parent;
	private Micrograph micrograph;
	
	
	
	public MicrographParticlesJDialog(XmippParticlePickerJFrame parent,
			Micrograph micrograph)
	{
		super(parent);
		this.parent = parent;
		this.micrograph = micrograph;
		initComponents();
	}
	
	private void initComponents()
	{
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Cells Info");
		JScrollPane sp = new JScrollPane();
		JList list = new JList();
		list.setOpaque(true);
		list.setModel(new ParticlesTableModel(micrograph.getParticles(), parent));
		//list.setPreferredSize(new Dimension(1000, 500));
		sp.setViewportView(list);
		add(sp);
		pack();
		centerScreen();
		setVisible(true);
	}
	private void centerScreen() {
		// TODO Auto-generated method stub
		
	}
	class ParticlesTableModel extends AbstractListModel
	{

		private List<Particle> particles;
		private XmippParticlePickerJFrame frame;

		public ParticlesTableModel(List<Particle> particles, XmippParticlePickerJFrame frame) {
			this.particles = particles;
			this.frame = frame;
		}

		
		@Override
		public Object getElementAt(int index) {
			Particle p = particles.get(index);
			return p.getImageIcon(micrograph.getImage(), p.getFamily().getSize());
		}

		@Override
		public int getSize() {
			return particles.size();
		}
		
	}

}
