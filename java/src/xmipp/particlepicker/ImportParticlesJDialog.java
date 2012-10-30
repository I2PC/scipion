package xmipp.particlepicker;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.io.File;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import xmipp.jni.Filename;
import xmipp.particlepicker.training.gui.TrainingPickerJFrame;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippWindowUtil;

public class ImportParticlesJDialog extends XmippDialog {

	ParticlePickerJFrame parent;
	private JTextField sourcetf;
	private JButton browsebt;
	private JComboBox jcbFormat;
	public Format format = Format.Auto;
	private XmippFileChooser xfc = null;
	protected String path;
	private String message;
	private boolean singleSelection = true;

	private static String[] FormatStrings = { "Automatic", "Xmipp 2.4",
			"Xmipp 3.0", "Eman" };
	private static Format[] FormatList = { Format.Auto, Format.Xmipp24,
			Format.Xmipp30, Format.Eman };

	public ImportParticlesJDialog(JFrame parent) {
		super(parent, "Import Particles", true);
		this.parent = (ParticlePickerJFrame) parent;
		xfc = new XmippFileChooser();
		xfc.setFileSelectionMode(XmippFileChooser.FILES_AND_DIRECTORIES);
		initComponents();
	}// constructor
	
	public void setMultiselectionEnabled(boolean value){
		singleSelection = !value;
		xfc.setMultiSelectionEnabled(value);
	}

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
		jcbFormat.addActionListener(this);
		panel.add(jcbFormat, XmippWindowUtil.getConstraints(gbc, 1, 0, 1));
		sourcetf = new JTextField(20);
		panel.add(sourcetf, XmippWindowUtil.getConstraints(gbc, 1, 1, 1));
		browsebt = XmippWindowUtil.getIconButton("folderopen.gif", this);
		panel.add(browsebt, XmippWindowUtil.getConstraints(gbc, 2, 1, 1));
	}// function createContent

	@Override
	public void handleActionPerformed(ActionEvent evt) {
		Object o = evt.getSource();
		if (o == browsebt)
			browseDirectory();
		else if (o == jcbFormat) {
			format = FormatList[jcbFormat.getSelectedIndex()];
		}
	}// function actionPerformed

	private void browseDirectory() {
		int returnVal = xfc.showOpenDialog(this);
		try {
			if (returnVal == XmippFileChooser.APPROVE_OPTION) {
				if (singleSelection)
					path = xfc.getSelectedPath();
				else {
					path = "";
					for (File f: xfc.getSelectedFiles())
						path += f.getPath() + " ";
				}
				path = path.trim();
				sourcetf.setText(path);
			}
		} catch (Exception ex) {
			showException(ex);
		}

	}

	@Override
	public void handleOk() {
		try {
			importParticles();
		} catch (Exception e) {
			XmippDialog.showException(parent, e);
		}
	}
	
	protected void importParticlesFromFile(){
		((TrainingPickerJFrame)parent).importParticlesFromFile(format, path);
	}
	
	private boolean existsSelectedPaths(){
		if (singleSelection)
			return Filename.exists(path);
		
		String[] parts = path.split(" ");		
		for (String s: parts)
			if (!Filename.exists(s))
				return false;
		return true;
	}//function existsSelectedPaths

	private void importParticles() {
		path = sourcetf.getText().trim();

		if (path == null || path.equals(""))
			showError(XmippMessage.getEmptyFieldMsg("directory"));
		else if (!existsSelectedPaths())
			showError(XmippMessage.getPathNotExistsMsg(path));
		else {			
			if (new File(path).isDirectory()) 
				parent.importParticlesFromFolder(format, path);
			else
				parent.importParticlesFromFile(format, path);
		}
	}
}
