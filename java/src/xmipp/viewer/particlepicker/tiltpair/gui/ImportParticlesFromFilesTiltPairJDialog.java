package xmipp.viewer.particlepicker.tiltpair.gui;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.NumberFormat;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import xmipp.utils.XmippMessage;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.ImportParticlesJDialog;


public class ImportParticlesFromFilesTiltPairJDialog extends ImportParticlesJDialog {

	

	private JTextField sourcetf2;
	private JButton browsebt2;
	private String path2;

	public ImportParticlesFromFilesTiltPairJDialog(TiltPairPickerJFrame parent) {
		super(parent);
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
		panel.add(new JLabel("Untilted Source:"),
				XmippWindowUtil.getConstraints(gbc, 0, 1, 1));
		
		panel.add(new JLabel("Tilted Source:"),
				XmippWindowUtil.getConstraints(gbc, 0, 2, 1));

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
		
		sourcetf2 = new JTextField(20);
		panel.add(sourcetf2, XmippWindowUtil.getConstraints(gbc, 1, 2, 1));
		browsebt = XmippWindowUtil.getIconButton("folderopen.gif", this);
		browsebt2 = XmippWindowUtil.getIconButton("folderopen.gif", this);
		browsebt2.setActionCommand("browse2");
		panel.add(browsebt, XmippWindowUtil.getConstraints(gbc, 2, 1, 1));
		panel.add(browsebt2, XmippWindowUtil.getConstraints(gbc, 2, 2, 1));
		panel.add(new JLabel("Scale To:"),
				XmippWindowUtil.getConstraints(gbc, 0, 3, 1));
		scaletf = new JFormattedTextField(NumberFormat.getNumberInstance());
		scaletf.setColumns(3);
		scaletf.setValue(1);
		panel.add(scaletf, XmippWindowUtil.getConstraints(gbc, 1, 3));
		
		panel.add(new JLabel("Invert X:"),
				XmippWindowUtil.getConstraints(gbc, 0, 4));
		invertxcb = new JCheckBox();
		panel.add(invertxcb, XmippWindowUtil.getConstraints(gbc, 1, 4));
		
		panel.add(new JLabel("Invert Y:"),
				XmippWindowUtil.getConstraints(gbc, 0, 5));
		invertycb = new JCheckBox();
		panel.add(invertycb, XmippWindowUtil.getConstraints(gbc, 1, 5));
		
		
	}// function createContent

	@Override
	public void handleActionPerformed(ActionEvent evt) {
		Object o = evt.getSource();
		if (o == browsebt)
			browseDirectory(sourcetf);
		else 
			browseDirectory(sourcetf2);
	}// function actionPerformed

	
	
	
	protected String importParticles() {
		path = sourcetf.getText().trim();

		if (path == null || path.equals("") )
			showError(XmippMessage.getEmptyFieldMsg("Untilted Source"));
		path2 = sourcetf2.getText().trim();

		if (path2 == null || path2.equals("") )
			showError(XmippMessage.getEmptyFieldMsg("Tilted Source"));
		return ((TiltPairPickerJFrame)parent).importParticlesFromFiles(format, path, path2, ((Number)scaletf.getValue()).floatValue(), invertxcb.isSelected(), invertycb.isSelected());
		
	}
	
}
