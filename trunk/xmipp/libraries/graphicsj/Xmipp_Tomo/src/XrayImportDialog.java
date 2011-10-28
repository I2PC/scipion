/***************************************************************************
 *
 * @author: Jesus Cuenca (jcuenca@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.util.Hashtable;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.io.SaveDialog;




/**
 * @author jcuenca
 *
 */

// TODO: refactor - reuse Form code for mainPanel
public class XrayImportDialog extends JDialog implements ActionListener
//, TextListener, FocusListener, ItemListener, KeyListener, AdjustmentListener, WindowListener 
{

	private JPanel mainPanel;
	private JPanel okCancelPanel;
	private JButton okButton,cancelButton;
	private GridBagConstraints gbc;
	//int mainPanelRow=0;
	private ExitValue status= ExitValue.CANCEL;
	private static String BROWSE_LABEL="Browse...", DATA_LABEL="Data directory ", FLAT_LABEL="Flatfield directory ", ROOT_LABEL="Root directory ", ROOT_FILE_LABEL="Root file ", CROP_LABEL="Crop ";
	private Hashtable<String, JTextField> textFields= new Hashtable<String, JTextField>();
	private String command;
	
	public XrayImportDialog(String title, Frame owner) {
		super(owner,title,true);
		setLayout(new BorderLayout());
		mainPanel=new JPanel();
		mainPanel.setLayout(new GridBagLayout());
		gbc=new GridBagConstraints();
		//gbc.ipadx=5; gbc.ipady=5;
		gbc.insets=new Insets(5,5,5,5);
		// gbc.gridwidth=3; gbc.gridheight=4;
		// mainPanel.setPreferredSize(new Dimension(3*100, 4*50));
		getContentPane().add(mainPanel,BorderLayout.CENTER);
		
		okCancelPanel=new JPanel();
		okCancelPanel.setLayout(new FlowLayout());
		okButton=new JButton("OK");
		cancelButton=new JButton("Cancel");
		okButton.addActionListener(this);
		cancelButton.addActionListener(this);
		okCancelPanel.add(okButton);
		okCancelPanel.add(cancelButton);
		getContentPane().add(okCancelPanel,BorderLayout.PAGE_END);
		
		
    	addBrowseField(0,DATA_LABEL, "/home/jcuenca/sample_data/20090606-10/10s/");
    	addBrowseField(1,FLAT_LABEL, "/home/jcuenca/sample_data/20090606-10/flatfields/");
    	addBrowseField(2,ROOT_LABEL, "/home/jcuenca/sample_data/20090606-10/ProcessedData");
    	addStringField(3,ROOT_FILE_LABEL, "img",10);
    	// addStringField(2,ROOT_LABEL, "ProcessedData/img",10);
    	// maybe it's better to use a slidebar instead of a numeric field?
    	addNumericField(4,CROP_LABEL, 7);
	}
	
      
    public void showDialog(){
    	pack();
    	setVisible(true);
    }
    
    // label identifies this field set
    // refactor like addTextField + browsefield specifics
    public void addBrowseField(int row, String label, String defaultText){

    	String label2 = label;
   		if (label2.indexOf('_')!=-1)
   			label2 = label2.replace('_', ' ');
		JLabel theLabel = new JLabel(label2);
		gbc.gridy = row; gbc.gridx=0;
		gbc.anchor = GridBagConstraints.EAST;
		mainPanel.add(theLabel,gbc);

		JTextField tf = new JTextField(defaultText, 15);
		tf.addActionListener(this);
		// tf.addTextListener(this);
		// tf.addFocusListener(this);
		// tf.addKeyListener(this);
		tf.setEditable(true);
		textFields.put(label, tf);
		gbc.gridx=1;
		mainPanel.add(tf,gbc);
		
		JButton b= new JButton(BROWSE_LABEL);
		b.setName(label);
		b.addActionListener(this);
		gbc.gridx=2;
		gbc.fill=GridBagConstraints.NONE;
		mainPanel.add(b,gbc);
    }
    
    /**
     * 
     * @param row (starts at 0)
     * @param label
     * @param defaultText
     * @param columns
     */
	public void addStringField(int row,String label, String defaultText, int columns) {
		JLabel theLabel = new JLabel(label);
		gbc.gridx = 0; gbc.gridy = row;
		gbc.weightx = 0.0; gbc.gridwidth=1;
		gbc.anchor = GridBagConstraints.EAST;

		mainPanel.add(theLabel,gbc);
		JTextField tf = new JTextField(defaultText, columns);
		textFields.put(label, tf);
		tf.addActionListener(this);
		gbc.gridx = 1; gbc.gridy = row;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.weightx = 0.0; gbc.gridwidth=GridBagConstraints.REMAINDER;
		tf.setEditable(true);
		mainPanel.add(tf,gbc);
    }
	
	public void addNumericField(int row,String label,int defaultValue){
		addStringField(row, label, String.valueOf(defaultValue), 3);
	}
	
	private void buildCommand(){
		setCommand("xmipp_xray_import --data " + getText(DATA_LABEL) + " --flat " + getText(FLAT_LABEL) + " --oroot " + getText(ROOT_LABEL)+
				"/" + getText(ROOT_FILE_LABEL) + " --crop " + getText(CROP_LABEL) + " --thr 1");
	}
	
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		String label = e.getActionCommand();
		if (source== okButton || source==cancelButton) {
			if(source == okButton)
				status=ExitValue.OK;
			buildCommand();
			dispose();
		}else if(BROWSE_LABEL.equals(label))
			actionBrowse(((JButton)source).getName());
		
		/* else
            notifyListeners(e);*/
	}


	private void actionBrowse(String fieldName){
		String path=TomoFileDialog.openDialog(fieldName, this);
		if(path != null) 
			textFields.get(fieldName).setText(path);
	}
	
	/**
	 * @return the status
	 */
	public ExitValue getStatus() {
		return status;
	}


	/**
	 * @param status the status to set
	 */
	private void setStatus(ExitValue status) {
		this.status = status;
	}

	

	/**
	 * Show a file open browser and...
	 * @deprecated - see FileDialog
	 * @return the path of the file chosen by the user, or empty string (not
	 *         null) if the user cancels the dialog
	 */
	public static String dialogOpen() {
		OpenDialog od = new OpenDialog("Import file", null);
		String directory = od.getDirectory();
		String fileName = od.getFileName();
		String path = "";
		if (fileName != null)
			path = directory + fileName;
		return path;
	}
	
	/**
	 * Show a file save browser and...
	 * 
	 * @return the path of the file chosen by the user, or empty string (not
	 *         null) if the user cancels the dialog
	 */
	public static String dialogSave(String initialPath,String extension) {
		SaveDialog sd = new SaveDialog("Save as...", initialPath,extension);
		String directory = sd.getDirectory();
		String fileName = sd.getFileName();
		String path = "";
		if (fileName != null)
			path = directory + fileName;
		return path;
	}
	
	/**
	 * @param title
	 *            Dialog title
	 * @param message
	 *            Dialog message
	 * @return Xmipp_Tomo.ExitValues.CANCEL / YES / NO
	 */
	
	// TODO: use swing JOptionPane dialog directly (no need of a static method...)
	public static ExitValue dialogYesNoCancel(String title, String message, boolean showCancelButton) {
		GenericDialog gd = new GenericDialog(title);

		gd.addMessage(message);
		if(showCancelButton)
			gd.enableYesNoCancel();
		gd.showDialog();
		if (gd.wasCanceled())
			return ExitValue.CANCEL;
		else if (gd.wasOKed())
			return ExitValue.YES;
		else
			return ExitValue.NO;
	}


	/**
	 * @return the command
	 */
	public  String getCommand() {
		return command;
	}
	
	/**
	 * 
	 * @return path to the mrcs result of the import
	 */
	public String getImagePath(){
		return getText(ROOT_LABEL) + "/" + getText(ROOT_FILE_LABEL) + ".mrcs";
	}


	/**
	 * @param command the command to set
	 */
	private void setCommand(String command) {
		this.command = command;
	}
	
	private String getText(String textField){
		String result ="";
		try{
			result= textFields.get(textField).getText();
		}catch (NullPointerException ex){
			result = "";
		}
		return result;
	}
	
	public static void main(String args[]) {
		JFrame window=new JFrame();
		Container content = window.getContentPane();
		//content.add(dialog, BorderLayout.CENTER);
		window.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {System.exit(0);}
		});

		window.setSize(200, 200);
		window.setVisible(true);
		XrayImportDialog dialog = new XrayImportDialog("Titulo", window);
		dialog.showDialog();

	}
}
