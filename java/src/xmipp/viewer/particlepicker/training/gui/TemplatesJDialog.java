package xmipp.viewer.particlepicker.training.gui;

import ij.ImagePlus;
import java.awt.Dimension;
import javax.swing.ImageIcon;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.ImageGeneric;
import xmipp.utils.XmippWindowUtil;

public class TemplatesJDialog extends JDialog
{

	protected SupervisedParticlePickerJFrame frame;
	protected JPanel templatespn;
	protected int width, height;

	public TemplatesJDialog(SupervisedParticlePickerJFrame frame)
	{
		super(frame);
		this.frame = frame;
		initComponents();
	}

	public synchronized void loadTemplates(boolean resize)
	{
		try
		{
			ImageGeneric templates = frame.getParticlePicker().getTemplates();
			int size = frame.getParticlePicker().getSize();
			templatespn.removeAll();
//			if (!frame.getParticlePicker().hasManualParticles())
				templatespn.setPreferredSize(new Dimension((int) (size * templates.getNDim()  + 20), size + 5));
			if(frame.getParticlePicker().hasManualParticles())
			{
				ImagePlus template;

				for (int i = 0; i < frame.getParticlePicker().getTemplatesNumber(); i++)
				{
					template = XmippImageConverter.convertToImagePlus(templates, ImageGeneric.FIRST_IMAGE + i);
					templatespn.add(new JLabel(new ImageIcon(template.getImage(), "")));

				}
			}
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();

		}
		templatespn.repaint();
		pack();
	}

	private void initComponents()
	{
		setDefaultCloseOperation(HIDE_ON_CLOSE);
		setTitle("Templates");
		templatespn = new JPanel();
		add(templatespn);
		loadTemplates(true);
		XmippWindowUtil.setLocation(0.6f, 0, this);
		setVisible(true);

	}

	public void close()
	{
		setVisible(false);
		dispose();

	}

}
