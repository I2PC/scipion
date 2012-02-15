package xmipp.particlepicker;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Panel;
import java.awt.ScrollPane;
import java.awt.event.ComponentEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.List;

import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import xmipp.particlepicker.tiltpair.gui.TiltPairParticlesJDialog;
import xmipp.particlepicker.tiltpair.model.TiltedParticle;
import xmipp.particlepicker.tiltpair.model.UntiltedParticle;
import xmipp.particlepicker.training.gui.TrainingPickerJFrame;
import xmipp.particlepicker.training.model.TrainingParticle;
import xmipp.utils.WindowUtils;
import xmipp.utils.XmippMessage;

public class ParticlesJDialog extends JDialog
{

	protected ParticlePickerJFrame frame;
	protected Panel particlespn;
	protected ScrollPane sp;
	protected GridBagConstraints constraints;
	protected int width, height;

	public ParticlesJDialog(ParticlePickerJFrame frame)
	{
		super(frame);
		this.frame = frame;
		initComponents();
		
		addWindowListener(new WindowAdapter()
		{
			public void windowClosing(WindowEvent winEvt)
			{
				resetParticlesJDialog();
			}

		});
	}

	protected void resetParticlesJDialog()
	{
		frame.particlesdialog = null;
		
	}

	public void loadParticles(boolean resize)
	{
		int side, rows, columns, width = 0, height = 0;
		List<? extends TrainingParticle> particles = frame.getParticles();
		side = frame.getSide(frame.getFamily().getSize());
		
		if(particles.isEmpty())
			throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("particles"));
		if(side == 0)
			throw new IllegalArgumentException(XmippMessage.getOutOfBoundsMsg("side"));
		
		if(resize)
		{
			columns = Math.min(200, particles.size() * side) / side;
			rows = (int) Math.ceil(particles.size() / (float) columns);
			width = side * columns;
			height = (side * Math.min(10, rows));
			boolean scroll = (height < rows * side);
			width = width + (scroll? 40: 20);
			height = height + (scroll? 0: 20);
			sp.setPreferredSize(new Dimension(width, height));
		}
		else
		{
			Dimension d = sp.getSize();
			columns = (int) d.getWidth() / side;
			rows = (int) Math.ceil((particles.size() / (float) columns));
		}
		particlespn.removeAll();
		particles = frame.getParticles();
		int index = 0;
		ParticleCanvas c;
		TrainingParticle p;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns; j++, index++)
			{
				if (index == particles.size())
					break;
				p = particles.get(index);
				c = p.getParticleCanvas(frame);
				particlespn.add(c, WindowUtils.getConstraints(constraints, j, i, 1));
			}
		if (resize)
			pack();
		particlespn.repaint();
	}

	private void initComponents()
	{
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Particles");
		constraints = new GridBagConstraints();
		sp = new ScrollPane();
		particlespn = new Panel(new GridBagLayout());
		sp.add(particlespn);
		add(sp);
		loadParticles(true);
		WindowUtils.setLocation(0.6, 0, this);
		setVisible(true);
		setAlwaysOnTop(true);
		this.addComponentListener(new java.awt.event.ComponentAdapter()
		{
			public void componentResized(ComponentEvent e)
			{
					loadParticles(false);
			}
		});
	}

	public void close()
	{
		setVisible(false);
		dispose();

	}

}
