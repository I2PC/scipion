package xmipp.viewer.particlepicker;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.text.NumberFormat;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import xmipp.jni.Filename;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.tiltpair.gui.ImportParticlesFromFilesTiltPairJDialog;
import xmipp.viewer.particlepicker.training.gui.SingleParticlePickerJFrame;

public class ImportParticlesJDialog extends XmippDialog {

	protected ParticlePickerJFrame parent;
	protected JTextField sourcetf;
	protected JButton browsebt;
	protected JComboBox jcbFormat;
	public Format format = Format.Auto;
	protected XmippFileChooser xfc = null;
	protected String path;
	private JFormattedTextField scaletf;
	private JCheckBox invertxcb;
	private JCheckBox invertycb;

	protected static String[] FormatStrings = { "Automatic", "Xmipp 2.4",
			"Xmipp 3.0", "Eman" };
	protected static Format[] FormatList = { Format.Auto, Format.Xmipp24,
			Format.Xmipp30, Format.Eman };

	public ImportParticlesJDialog(ParticlePickerJFrame parent) {
		super(parent, "Import Particles", true);
		this.parent = parent;
		xfc = new XmippFileChooser();
		if(parent instanceof SingleParticlePickerJFrame)
			xfc.setFileSelectionMode(XmippFileChooser.FILES_AND_DIRECTORIES);
		else if (this instanceof ImportParticlesFromFilesTiltPairJDialog)
			xfc.setFileSelectionMode(XmippFileChooser.FILES_ONLY);
		else 
			xfc.setFileSelectionMode(XmippFileChooser.DIRECTORIES_ONLY);
		
		xfc.setMultiSelectionEnabled(false);
		initComponents();
	}// constructor
	
	

	@Override
	protected void createContent(JPanel panel) {
		setResizable(false);
		panel.setLayout(new GridBagLayout());
		GridBagConstraints gbc = new GridBagConstraints();
		gbc.insets = new Insets(5, 5, 5, 5);

		gbc.anchor = GridBagConstraints.EAST;
		panel.add(new JLabel("Format:"),
				XmippWindowUtil.getConstraints(gbc, 0, 0, 1));
		panel.add(new JLabel("Source:"),
				XmippWindowUtil.getConstraints(gbc, 0, 1, 1));

		gbc.anchor = GridBagConstraints.WEST;
		/** Create a combobox with possible formats */
		jcbFormat = new JComboBox(FormatStrings);
		jcbFormat.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent arg0) {
				format = FormatList[jcbFormat.getSelectedIndex()];
				
			}
		});
		panel.add(jcbFormat, XmippWindowUtil.getConstraints(gbc, 1, 0, 1));
		
		sourcetf = new JTextField(20);
		panel.add(sourcetf, XmippWindowUtil.getConstraints(gbc, 1, 1, 1));
		
		browsebt = XmippWindowUtil.getIconButton("folderopen.gif", this);
		panel.add(browsebt, XmippWindowUtil.getConstraints(gbc, 2, 1, 1));
		
		panel.add(new JLabel("Scale To:"),
				XmippWindowUtil.getConstraints(gbc, 0, 2, 1));
		scaletf = new JFormattedTextField(NumberFormat.getNumberInstance());
		scaletf.setColumns(3);
		scaletf.setValue(1);
		panel.add(scaletf, XmippWindowUtil.getConstraints(gbc, 1, 2));
		
		panel.add(new JLabel("Invert X:"),
				XmippWindowUtil.getConstraints(gbc, 0, 3));
		invertxcb = new JCheckBox();
		panel.add(invertxcb, XmippWindowUtil.getConstraints(gbc, 1, 3));
		
		panel.add(new JLabel("Invert Y:"),
				XmippWindowUtil.getConstraints(gbc, 0, 4));
		invertycb = new JCheckBox();
		panel.add(invertycb, XmippWindowUtil.getConstraints(gbc, 1, 4));
		
	}// function createContent

	@Override
	public void handleActionPerformed(ActionEvent evt) {
			browseDirectory(sourcetf);
		
	}// function actionPerformed

	protected void browseDirectory(JTextField sourcetf) {
		int returnVal = xfc.showOpenDialog(this);
		try {
			if (returnVal == XmippFileChooser.APPROVE_OPTION) {
				path = xfc.getSelectedPath();
				sourcetf.setText(path);
			}
		} catch (Exception ex) {
			showException(ex);
		}

	}

	@Override
	public void handleOk() {
		try {
			path = sourcetf.getText().trim();

			if (path == null || path.equals(""))
				showError(XmippMessage.getEmptyFieldMsg("source"));
			else if (!existsSelectedPath())
				showError(XmippMessage.getPathNotExistsMsg(path));
			else
			{
				String result = parent.importParticles(format, path, ((Number)scaletf.getValue()).floatValue(), invertxcb.isSelected(), invertycb.isSelected());
				if(result != null && !result.isEmpty())
					JOptionPane.showMessageDialog(this, result);
			}
		} catch (Exception e) {
			XmippDialog.showException(parent, e);
		}
	}
	
	
	
	
	private boolean existsSelectedPath(){
			return Filename.exists(path);
		
	}//function existsSelectedPaths

	
	


}
