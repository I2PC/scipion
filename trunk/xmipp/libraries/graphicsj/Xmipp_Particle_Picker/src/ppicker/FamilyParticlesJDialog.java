package ppicker;

import ij.ImagePlus;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.LayoutManager;
import java.util.List;

import javax.swing.AbstractListModel;
import javax.swing.JDialog;
import javax.swing.JList;
import javax.swing.JScrollPane;


public class FamilyParticlesJDialog extends JDialog {
	
	private XmippParticlePickerJFrame parent;
	private List<Particle> particles;
	
	
	
	public FamilyParticlesJDialog(XmippParticlePickerJFrame parent, List<Particle> particles)
	{
		super(parent);
		this.parent = parent;
		this.particles = particles;
		initComponents();
	}
	
	private void initComponents()
	{
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Cells Info");
		JScrollPane sp = new JScrollPane();
		JList list = new JList();
		list.setOpaque(true);
		list.setModel(new ParticlesTableModel(particles, parent));
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
			return p.getImageIcon(frame.getImage(), p.getFamily().getRadius());
		}

		@Override
		public int getSize() {
			return particles.size();
		}
		
	}

}
