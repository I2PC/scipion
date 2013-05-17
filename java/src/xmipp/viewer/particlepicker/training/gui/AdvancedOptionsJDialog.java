package xmipp.viewer.particlepicker.training.gui;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.ImageCanvas;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.Insets;
import java.awt.Panel;
import java.awt.ScrollPane;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.text.NumberFormat;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import sun.net.www.content.image.jpeg;

import xmipp.ij.commons.ImagePlusLoader;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.jni.ImageGeneric;
import xmipp.utils.XmippWindowUtil;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;

public class AdvancedOptionsJDialog extends JDialog {

	protected SingleParticlePickerJFrame frame;
	protected JPanel templatespn;
	protected int width, height;
	private JFormattedTextField templatestf;
	private JButton loadtemplatesbt;
	private JLabel checkpercentlb;
	private JFormattedTextField autopickpercenttf;
	private JButton okbt;

	public AdvancedOptionsJDialog(SingleParticlePickerJFrame frame) {
		super(frame);
		this.frame = frame;

		if(!frame.getParticlePicker().hasManualParticles())
			throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("Particles"));

		initComponents();

		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent winEvt) {
				resetTemplatesJDialog();
			}

		});
	}

	protected void resetTemplatesJDialog() {
		frame.optionsdialog = null;

	}

	public void loadTemplates(boolean resize) {

		try {
			ImageGeneric templates = frame.getParticlePicker().getTemplates();
			
			int size = frame.getParticlePicker().getSize();

			if (!frame.getParticlePicker().hasManualParticles()) {

				templatespn.removeAll();
				templatespn.setPreferredSize(new Dimension(
						(int) (size * templates.getNDim()), size));
				pack();
				return;
			}

			templatespn.removeAll();
			ImagePlus template;
			for (int i = 0; i < frame.getParticlePicker().getTemplatesNumber(); i ++) {
				template = frame.getParticlePicker().getTemplatesImage(ImageGeneric.FIRST_IMAGE + i);
				templatespn.add(new ImageCanvas(template));

			}
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		templatespn.repaint();
		pack();
	}

	private void initComponents() {
		setDefaultCloseOperation(HIDE_ON_CLOSE);
		setTitle("Advanced Options");
		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(0, 5, 0, 5);
		constraints.anchor = GridBagConstraints.WEST;
		setLayout(new GridBagLayout());
		
		add(new JLabel("Templates:"), XmippWindowUtil.getConstraints(constraints, 0, 0));

		templatestf = new JFormattedTextField(NumberFormat.getNumberInstance());
		templatestf.setColumns(3);
		templatestf.setValue(frame.getParticlePicker().getTemplatesNumber());
		loadtemplatesbt = XmippWindowUtil.getTextButton("View", new ActionListener()
		{
			
			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				loadTemplates();
			}
		});
		
		add(templatestf, XmippWindowUtil.getConstraints(constraints, 1, 0));
		add(templatestf);
		templatespn = new JPanel();
		add(loadtemplatesbt, XmippWindowUtil.getConstraints(constraints, 2, 0));

		checkpercentlb = new JLabel("Autopick Check (%):");
		add(checkpercentlb, XmippWindowUtil.getConstraints(constraints, 0, 1));
		autopickpercenttf = new JFormattedTextField(NumberFormat.getIntegerInstance());
		autopickpercenttf.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				if (autopickpercenttf.getValue() == null)
				{
					JOptionPane.showMessageDialog(AdvancedOptionsJDialog.this, XmippMessage.getEmptyFieldMsg("Check (%)"));
					autopickpercenttf.setValue(frame.getMicrograph().getAutopickpercent());
					return;
				}

				int autopickpercent = ((Number) autopickpercenttf.getValue()).intValue();
				frame.getMicrograph().setAutopickpercent(autopickpercent);
				frame.getParticlePicker().setAutopickpercent(autopickpercent);
				

			}
		});
		autopickpercenttf.setColumns(3);
		add(autopickpercenttf, XmippWindowUtil.getConstraints(constraints, 1, 1));
		
		okbt = XmippWindowUtil.getTextButton("Ok", new ActionListener()
		{
			
			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				setVisible(false);
				dispose();
				
			}
		});
		add(okbt, XmippWindowUtil.getConstraints(constraints, 2, 2));
		loadTemplates(true);
		XmippWindowUtil.setLocation(0.6f, 0, this);
		setVisible(true);
		setAlwaysOnTop(true);
		// this.addComponentListener(new java.awt.event.ComponentAdapter() {
		// public void componentResized(ComponentEvent e) {
		// loadTemplates(false);
		// }
		// });
	}

	protected void loadTemplates()
	{
		// TODO Auto-generated method stub
		
	}

	public void close() {
		setVisible(false);
		dispose();

	}

}
