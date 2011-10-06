package trainingpicker.gui;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Panel;
import java.awt.ScrollPane;
import java.awt.event.AdjustmentEvent;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowStateListener;
import java.util.List;

import javax.swing.JDialog;

import trainingpicker.model.TrainingParticle;

public class ParticlesJDialog extends JDialog
{

	private ParticlePickerJFrame frame;
	private List<TrainingParticle> particles;
	private int rows, columns;
	private int side;
	private Panel particlespn;
	private ScrollPane sp;
	private GridBagConstraints constraints;
	private int width;
	private int height;

	public ParticlesJDialog(ParticlePickerJFrame frame)
	{
		super(frame);
		this.frame = frame;
		initComponents();
		
	}
	
	public void loadParticles()
	{
		side = (int)(frame.getFamily().getSize() * frame.getMagnification());
		columns = Math.min(800, frame.getFamilyData().getParticles().size() * side) / side;
		rows = (int)Math.ceil((frame.getFamilyData().getParticles().size()/(float)columns));
		width = side * columns;
		height = (side * Math.min(3, rows));
		particlespn.removeAll();
		particles = frame.getFamilyData().getParticles();
		int index = 0;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns; j++, index++)
			{
				if (index == particles.size())
					break;
				particlespn.add(particles.get(index).getParticleCanvas(frame), WindowUtils.getConstraints(constraints, j, i, 1));
			}
		sp.setPreferredSize(new Dimension(width + 20, height));
		pack();
	}

	private void initComponents()
	{
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Particles");
		setResizable(false);
		constraints = new GridBagConstraints();
		sp = new ScrollPane();
		particlespn = new Panel(new GridBagLayout());
		sp.add(particlespn);
		add(sp);
		
		loadParticles();
//		addComponentListener(new ComponentListener()
//		{
//			
//			@Override
//			public void componentShown(ComponentEvent e)
//			{
//				// TODO Auto-generated method stub
//				
//			}
//			
//			@Override
//			public void componentResized(ComponentEvent e)
//			{
//				
//				Dimension d = sp.getSize();
//				columns = (int)d.getWidth()/side;
//				rows = (int)Math.ceil((particles.size()/(float)columns));
//				width = (int)d.getWidth();
//				height = (int)d.getHeight();
//				loadParticles();
//			}
//			
//			@Override
//			public void componentMoved(ComponentEvent e)
//			{
//				// TODO Auto-generated method stub
//				
//			}
//			
//			@Override
//			public void componentHidden(ComponentEvent e)
//			{
//				// TODO Auto-generated method stub
//				
//			}
//		});
		addWindowListener(new WindowAdapter()
		{
			public void windowClosing(WindowEvent winEvt)
		    {
				frame.particlesdialog = null;
		    }
		});
		WindowUtils.centerScreen(0, 0.9, this);
		setVisible(true);
	}

	public void close()
	{
		setVisible(false);
		dispose();
		
	}

	

}
