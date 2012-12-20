package xmipp.particlepicker.training.gui;

import ij.ImagePlus;
import ij.gui.ImageCanvas;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.Panel;
import java.awt.ScrollPane;
import java.awt.event.ComponentEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JDialog;
import javax.swing.JLabel;

import xmipp.ij.commons.ImagePlusLoader;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.jni.ImageGeneric;
import xmipp.particlepicker.ParticlePickerJFrame;
import xmipp.particlepicker.training.model.TrainingParticle;
import xmipp.utils.XmippWindowUtil;
import xmipp.utils.XmippMessage;

public class TemplatesJDialog extends JDialog {

	protected TrainingPickerJFrame frame;
	protected Panel templatespn;
	protected int width, height;

	public TemplatesJDialog(TrainingPickerJFrame frame) {
		super(frame);
		this.frame = frame;
		initComponents();

		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent winEvt) {
				resetTemplatesJDialog();
			}

		});
	}

	protected void resetTemplatesJDialog() {
		frame.templatesdialog = null;

	}

	public void loadTemplates(boolean resize) {

		try {
			frame.updateTemplates();
			ImageGeneric templates = frame.getFamily().getTemplates();
			
			int size = frame.getFamily().getSize();


			if (!frame.getParticlePicker().hasParticles()) {
				templatespn.removeAll();
				templatespn.setPreferredSize(new Dimension(
						(int) (size * templates.getNDim()), size));
				pack();
				return;
			}

			templatespn.removeAll();
			ImagePlus template;

			for (int index = 0; index < templates.getNDim(); index++) {
				template = XmippImageConverter.convertToImagePlus(templates,
						ImageGeneric.FIRST_IMAGE + index);
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
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Templates");
		templatespn = new Panel();
		add(templatespn);
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

	public void close() {
		setVisible(false);
		dispose();

	}

}
