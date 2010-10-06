import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.Panel;
import java.awt.TextField;
import java.util.Hashtable;
import java.util.Vector;
import java.awt.Button;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentListener;
import java.awt.event.FocusListener;
import java.awt.event.ItemListener;
import java.awt.event.KeyListener;
import java.awt.event.TextListener;
import java.awt.event.WindowListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.io.SaveDialog;
import ij.plugin.frame.Recorder;



/**
 * @author jcuenca
 *
 */
public class XrayImportDialog extends JDialog implements ActionListener
//, TextListener, FocusListener, ItemListener, KeyListener, AdjustmentListener, WindowListener 
{

	private JPanel mainPanel;
	private JPanel okCancelPanel;
	private JButton okButton,cancelButton;
	private GridBagConstraints gbc;
	//int mainPanelRow=0;
	private Xmipp_Tomo.ExitValues status= Xmipp_Tomo.ExitValues.CANCEL;
	private static String BROWSE_LABEL="Browse...", DATA_LABEL="Data", FLAT_LABEL="Flat", ROOT_LABEL="Root", CROP_LABEL="Crop";
	private Hashtable<String, JTextField> textFields= new Hashtable<String, JTextField>();
	private String command;
	
	public XrayImportDialog(String title, Frame owner) {
		super(owner,title,true);
		// TODO Auto-generated constructor stub
		setLayout(new BorderLayout());
		mainPanel=new JPanel();
		mainPanel.setLayout(new GridBagLayout());
		gbc=new GridBagConstraints();
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
		
		
    	addBrowseField(0,DATA_LABEL, "");
    	addBrowseField(1,FLAT_LABEL, "");
    	addStringField(2,ROOT_LABEL, "ProcessedData/img",10);
    	// maybe it's better to use a slidebar instead of a numeric field?
    	addNumericField(3,CROP_LABEL, 7);
	}
	
    
    protected JLabel makeLabel(String label) {
    	if (IJ.isMacintosh())
    		label += " ";
		return new JLabel(label);
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
		JLabel theLabel = makeLabel(label2);
		gbc.gridy = row; gbc.gridx=0;
		mainPanel.add(theLabel,gbc);

		JTextField tf = new JTextField(defaultText, 15);
		if (IJ.isLinux()) tf.setBackground(Color.white);
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
		mainPanel.add(b,gbc);
		//mainPanelRow++;
    }
    
	public void addStringField(int row,String label, String defaultText, int columns) {
   		String label2 = label;
   		if (label2.indexOf('_')!=-1)
   			label2 = label2.replace('_', ' ');
		JLabel theLabel = makeLabel(label2);
		gbc.gridx = 0; gbc.gridy = row;
		gbc.anchor = GridBagConstraints.EAST;
		//c.gridwidth = 1;
		/*boolean custom = customInsets;
		if (stringField==null) {
			stringField = new Vector(4);
			c.insets = getInsets(5, 0, 5, 0);
		} else
			c.insets = getInsets(0, 0, 5, 0);
		grid.setConstraints(theLabel, c);*/
		mainPanel.add(theLabel,gbc);
		/*if (custom) {
			if (stringField.size()==0)
				c.insets = getInsets(5, 0, 5, 0);
			else
				c.insets = getInsets(0, 0, 5, 0);
		}*/
		JTextField tf = new JTextField(defaultText, columns);
		if (IJ.isLinux()) tf.setBackground(Color.white);
		/*tf.setEchoChar(echoChar);
		echoChar = 0;*/
		textFields.put(label, tf);
		tf.addActionListener(this);
		/*tf.addTextListener(this);
		tf.addFocusListener(this);
		tf.addKeyListener(this);*/
		gbc.gridx = 1; gbc.gridy = row;
		gbc.anchor = GridBagConstraints.WEST;
		//grid.setConstraints(tf, c);
		tf.setEditable(true);
		mainPanel.add(tf,gbc);
		/*stringField.addElement(tf);
		if (Recorder.record || macro)
			saveLabel(tf, label);
		y++;*/
    }
	
	public void addNumericField(int row,String label,int defaultValue){
		addStringField(row, label, String.valueOf(defaultValue), 3);
	}
	
	private void buildCommand(){
		setCommand("xmipp_ray -import -data " + getText(DATA_LABEL) + " -flat " + getText(FLAT_LABEL) + "-oroot " + getText(ROOT_LABEL)+
				"-crop " + getText(CROP_LABEL) + "-thr 1");
	}
	
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		String label = e.getActionCommand();
		if (source== okButton || source==cancelButton) {
			if(source == okButton)
				status=Xmipp_Tomo.ExitValues.OK;
			buildCommand();
			dispose();
		}else if(BROWSE_LABEL.equals(label))
			actionBrowse(((JButton)source).getName());
		
		/* else
            notifyListeners(e);*/
	}


	private void actionBrowse(String fieldName){
		String path=dialogOpen();
		textFields.get(fieldName).setText(path);
	}
	
	/**
	 * @return the status
	 */
	public Xmipp_Tomo.ExitValues getStatus() {
		return status;
	}


	/**
	 * @param status the status to set
	 */
	private void setStatus(Xmipp_Tomo.ExitValues status) {
		this.status = status;
	}

	

	/**
	 * Show a file open browser and...
	 * 
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
	public static Xmipp_Tomo.ExitValues dialogYesNoCancel(String title, String message) {
		GenericDialog gd = new GenericDialog(title);

		gd.addMessage(message);
		gd.enableYesNoCancel();
		gd.showDialog();
		if (gd.wasCanceled())
			return Xmipp_Tomo.ExitValues.CANCEL;
		else if (gd.wasOKed())
			return Xmipp_Tomo.ExitValues.YES;
		else
			return Xmipp_Tomo.ExitValues.NO;
	}


	/**
	 * @return the command
	 */
	public  String getCommand() {
		return command;
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
}
