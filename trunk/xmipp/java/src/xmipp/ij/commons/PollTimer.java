package xmipp.ij.commons;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Timer;

public class PollTimer extends Timer
{
	private final XmippIJWindow xp;
	private static final int PERIOD = 5000;  // repeat every 5 seconds

	public PollTimer(final XmippIJWindow xp)
	{
		super(PERIOD, new ActionListener()
		{
			
			@Override
			public void actionPerformed(ActionEvent e)
			{
				xp.loadData();
			}
		});
		this.xp = xp;
	}
}
