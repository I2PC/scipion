package particlepicker.tiltpair.gui;

import java.awt.Dimension;
import java.util.List;

import particlepicker.ParticlePickerJFrame;
import particlepicker.ParticlesJDialog;
import particlepicker.WindowUtils;
import particlepicker.tiltpair.model.UntiltedParticle;
import particlepicker.training.gui.ParticleCanvas;
import particlepicker.training.model.Constants;
import particlepicker.training.model.TrainingParticle;

public class TiltPairParticlesJDialog extends ParticlesJDialog
{

	public TiltPairParticlesJDialog(ParticlePickerJFrame frame)
	{

		super(frame);

	}

	public void loadParticles(boolean resize)
	{
		int side, rows, columns, width = 0, height = 0;
		List<? extends TrainingParticle> particles = frame.getParticles();
		side = (int) (frame.getFamily().getSize() * frame.getMagnification());

		if (particles.isEmpty())
			throw new IllegalArgumentException(Constants.getEmptyFieldMsg("particles"));
		if (side == 0)
			throw new IllegalArgumentException(Constants.getOutOfBoundsMsg("side"));

		if (resize)
		{
			columns = Math.min(200, particles.size() * side * 2) / (side * 2);
			rows = (int) Math.ceil(particles.size() / (float) columns);
			width = side * columns * 2;
			height = (side * Math.min(10, rows));
			boolean scroll = (height < rows * side);
			width = width + (scroll ? 40 : 20);
			height = height + (scroll ? 0 : 20);
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
		UntiltedParticle up;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns * 2; j+= 2, index ++)
			{
				if (index == particles.size())
					break;
				p = particles.get(index);
				c = p.getParticleCanvas(frame);
				up = (UntiltedParticle) p;
				particlespn.add(c, WindowUtils.getConstraints(constraints, j, i, 1));
				if (up.getTiltedParticle() != null)
					particlespn.add(up.getTiltedParticle().getParticleCanvas(frame), WindowUtils.getConstraints(constraints, j + 1, i, 1));
			}
		if (resize)
			pack();
	}
}
