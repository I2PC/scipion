package xmipp.viewer.particlepicker;

import java.awt.BorderLayout;
import java.awt.Dialog;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Panel;
import java.awt.ScrollPane;
import java.awt.event.ComponentEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.List;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippWindowUtil;

public class ParticlesDialog extends Dialog
{

	protected ParticlePickerJFrame frame;
	protected Panel particlespn;
	protected ScrollPane sp;
	protected GridBagConstraints constraints;
	protected int width, height, rows, columns, side;
	public static final int maxHeight = (int)(0.9 * XmippWindowUtil.getScreenRectangle().height);
	public static final int minWidth = 600;

	public ParticlesDialog(ParticlePickerJFrame frame)
	{
		super(frame);
		this.frame = frame;
		initComponents();

		addWindowListener(new WindowAdapter()
		{
			public void windowClosing(WindowEvent winEvt)
			{
				resetParticlesDialog();
                close();
			}

		});
	}

	protected void resetParticlesDialog()
	{
		frame.particlesdialog = null;

	}

	public void loadParticles(boolean resize)
	{

		List<? extends PickerParticle> particles = frame.getAvailableParticles();
		side = frame.getSide();

		if (side == 0)
			throw new IllegalArgumentException(XmippMessage.getOutOfBoundsMsg("side"));

		if (particles.isEmpty())
		{
			particlespn.removeAll();
            width = minWidth;
            height = maxHeight;
			sp.setPreferredSize(new Dimension(width, height));
			pack();
			return;
		}

		if (resize) // first time or keep size
		{
			columns = Math.min(minWidth, particles.size() * side) / side;
			rows = (int) Math.ceil(particles.size() / (float) columns);
			width = side * columns + 20;
			height = Math.min(side * rows, maxHeight);
			
		}
		else//size changed by user
		{
			Dimension d = sp.getSize();
            width = (int) d.getWidth();
            height = (int)d.getHeight();
			columns = width / side;
			rows = (int) Math.ceil((particles.size() / (float) columns));
                        
		}
        sp.setPreferredSize(new Dimension(width, height));
		particlespn.removeAll();
                
		particles = frame.getAvailableParticles();
		int index = 0;
		ParticleCanvas c;
		PickerParticle p;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns; j++, index++)
			{
                                
				if (index == particles.size())
					break;
				p = particles.get(index);
				c = p.getParticleCanvas(frame);
                                
				particlespn.add(c, XmippWindowUtil.getConstraints(constraints, j, i, 1));
			}
        // particlespn.revalidate();
        sp.getVAdjustable().setValue(sp.getVAdjustable().getMaximum());
        sp.getHAdjustable().setValue(sp.getHAdjustable().getMaximum());
		pack();

	}
        
        

	private void initComponents()
	{
		
		setTitle("Particles");
        setLayout(new BorderLayout());
		constraints = new GridBagConstraints();
		sp = new ScrollPane(ScrollPane.SCROLLBARS_AS_NEEDED);
                
		particlespn = new Panel(new GridBagLayout());
		sp.add(particlespn);
		add(sp, BorderLayout.CENTER);
		loadParticles(true);
		XmippWindowUtil.setLocation(0.6f, 0, this);
		setVisible(true);
		setAlwaysOnTop(true);
		sp.addComponentListener(new java.awt.event.ComponentAdapter()
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
