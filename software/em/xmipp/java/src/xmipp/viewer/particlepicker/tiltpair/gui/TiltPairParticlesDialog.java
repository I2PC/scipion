package xmipp.viewer.particlepicker.tiltpair.gui;

import java.awt.Dimension;
import java.util.List;

import xmipp.utils.XmippMessage;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.ParticleCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.ParticlesDialog;
import xmipp.viewer.particlepicker.PickerParticle;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedParticle;

public class TiltPairParticlesDialog extends ParticlesDialog
{
    

	public TiltPairParticlesDialog(ParticlePickerJFrame frame)
	{

		super(frame);

	}

	public void loadParticles(boolean changesize)
	{
		
		List<? extends PickerParticle> particles = frame.getAvailableParticles();
		side = frame.getSide();

		if (side == 0)
			throw new IllegalArgumentException(XmippMessage.getOutOfBoundsMsg("side"));
		
		if (particles.isEmpty())
		{
			particlespn.removeAll();
			width = side * 6 + 10;
            height = maxHeight;
			sp.setPreferredSize(new Dimension(width, height));
			pack();
			return;
		}

		if (changesize)// first time or keep size
		{
			columns = 3;//Math.min(200, particles.size() * side * 2) / (side * 2);
			rows = (int) Math.ceil(particles.size() / (float) columns);
			width = side * columns * 2;
			height = Math.min(maxHeight, side * rows);
			boolean scroll = (height < rows * side);
			width = width + (scroll ? 40 : 20);
			height = height + (scroll ? 0 : 20);
			sp.setPreferredSize(new Dimension(width, height));
		}
		else
		{
			Dimension d = sp.getSize();
            width = (int) d.getWidth();
            height = (int)d.getHeight();
			columns = (int) d.getWidth() / (side * 2);
			rows = (int) Math.ceil((particles.size() / (float) columns));
		}
                sp.setPreferredSize(new Dimension(width, height));
                
		particlespn.removeAll();
		particles = frame.getAvailableParticles();
		int index = 0;
		ParticleCanvas c;
		PickerParticle p;
		UntiltedParticle up;
		
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns * 2; j+= 2, index ++)
			{
				if (index == particles.size())
					break;
				p = particles.get(index);
				c = p.getParticleCanvas(frame);
				up = (UntiltedParticle) p;
				particlespn.add(c, XmippWindowUtil.getConstraints(constraints, j, i, 1));
				
				if (up.getTiltedParticle() != null)
					particlespn.add(up.getTiltedParticle().getParticleCanvas(frame), XmippWindowUtil.getConstraints(constraints, j + 1, i, 1));
			}
		
                sp.getVAdjustable().setValue(sp.getVAdjustable().getMaximum());
                sp.getHAdjustable().setValue(sp.getHAdjustable().getMaximum());
                
		 pack();
		
	}
}
