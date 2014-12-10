package xmipp.viewer.particlepicker.training.gui;

import ij.ImagePlus;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.NumberFormat;
import javax.swing.ImageIcon;
import javax.swing.JDialog;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.ImageGeneric;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.training.model.Mode;

public class TemplatesJDialog extends JDialog
{

	protected SupervisedPickerJFrame frame;
	protected JPanel templatespn;
	protected int width, height;
    private JFormattedTextField templatestf;

	public TemplatesJDialog(SupervisedPickerJFrame frame)
	{
		super(frame);
		this.frame = frame;
		initComponents();
	}

	public void loadTemplates(boolean resize)
	{
		try
		{
			ImageGeneric templates = frame.getParticlePicker().getTemplates();
			int size = frame.getParticlePicker().getSize();
			templatespn.removeAll();
			templatespn.setPreferredSize(new Dimension((int) (size * templates.getNDim()  + 50), size + 5));

			if(frame.getParticlePicker().hasManualParticles())
			{
				ImagePlus template;
				for (int i = 0; i < frame.getParticlePicker().getTemplatesNumber(); i ++)
				{
					template = XmippImageConverter.convertToImagePlus(templates, ImageGeneric.FIRST_IMAGE + i);
					templatespn.add(new JLabel(new ImageIcon(template.getImage()), SwingConstants.CENTER));
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
                GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(0, 5, 0, 5);
		constraints.anchor = GridBagConstraints.WEST;
		setLayout(new GridBagLayout());
                add(new JLabel("Templates:"), XmippWindowUtil.getConstraints(constraints, 0, 0));

		templatestf = new JFormattedTextField(NumberFormat.getNumberInstance());
		templatestf.setColumns(3);
		templatestf.setValue(frame.getParticlePicker().getTemplatesNumber());

		templatestf.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				setTemplates();
			}
		});
                templatestf.setEnabled(frame.getParticlePicker().getMode() == Mode.Manual);
		add(templatestf, XmippWindowUtil.getConstraints(constraints, 1, 0));
                templatespn = new JPanel();
		add(templatespn, XmippWindowUtil.getConstraints(constraints, 0, 1, 2, GridBagConstraints.HORIZONTAL));
                
		loadTemplates(true);

		XmippWindowUtil.setLocation(0.6f, 0, this);
		setVisible(true);

	}

	public void close()
	{
		setVisible(false);
		dispose();

	}
        
        protected void setTemplates()
	{
		if (templatestf.getValue() == null || ((Number)templatestf.getValue()).intValue() <= 0)
		{
			XmippDialog.showInfo(frame, XmippMessage.getOutOfBoundsMsg("Templates"));
			templatestf.setValue(frame.getParticlePicker().getTemplatesNumber());
			return;
		}

		int templates = ((Number) templatestf.getValue()).intValue();
		if (templates != frame.getParticlePicker().getTemplatesNumber())
		{
			frame.getParticlePicker().setTemplatesNumber(templates);
			frame.loadTemplates();
		}
		
	}


}
