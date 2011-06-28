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
/**
 * Why? To encapsulate the peculiarities of GridBagLayout,
 * simplifying the creating of standard forms
 */
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.util.Hashtable;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

public class Form extends JPanel 
{

	private JPanel mainPanel;
	private GridBagConstraints gbc;
	private Hashtable<String, JTextField> textFields= new Hashtable<String, JTextField>();

	private static String BROWSE_LABEL="Browse...";

	public Form(String title) {
		super();
		setLayout(new BorderLayout());
		mainPanel=new JPanel();
		mainPanel.setLayout(new GridBagLayout());
		gbc=new GridBagConstraints();
		gbc.insets=new Insets(5,5,5,5);
		add(mainPanel,BorderLayout.CENTER);
	}

	public void addStringField(int row,String label, String defaultText) {
		addStringField(row,label, defaultText,true);
	}
	
	// TODO: bug? some labels don't align to right
	private void addStringField(int row,String label, String defaultText,boolean expandTextField) {
		JLabel theLabel = new JLabel(label);
		gbc.gridx = 0; gbc.gridy = row;
		gbc.weightx = 0.0; gbc.gridwidth=1;
		gbc.anchor = GridBagConstraints.EAST;
		mainPanel.add(theLabel,gbc);

		JTextField tf = new JTextField(defaultText);
		textFields.put(label, tf);
		gbc.gridx = 1; gbc.gridy = row;
		gbc.anchor = GridBagConstraints.WEST;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.weightx = 1.0;
		if(expandTextField)
			gbc.gridwidth=GridBagConstraints.REMAINDER;		
		tf.setEditable(true);
		mainPanel.add(tf,gbc);
	}

	public void addNumericField(int row,String label,int defaultValue){
		addStringField(row, label, String.valueOf(defaultValue));
	}

	public void addBrowseField(int row, String label, String defaultText){
		addStringField(row,label,defaultText,false);
		Command browseCommand = new Command("form.browse"+String.valueOf(row),BROWSE_LABEL, "openDialog", true, null);		
		Action action= new Action(this, browseCommand,label);
		JButton b = new JButton(action);
		gbc.gridx=2;
		gbc.fill=GridBagConstraints.NONE;
		mainPanel.add(b,gbc);
	}

	public void openDialog(String fieldName){
		String path=FileDialog.openDialog(fieldName, this);
		if(path != null) 
			textFields.get(fieldName).setText(path);
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
		JFrame window=new JFrame("Form test");
		Container content = window.getContentPane();
		content.setForeground(Color.BLUE);
		Form form = new Form("Title");
		form.setBackground(Color.GREEN);
		form.addStringField(0, "First field", "default");
		form.addStringField(1, "Second field", "");
		form.addBrowseField(2, "File", "/home/abc/file.txt");
		content.add(form, BorderLayout.CENTER);
		window.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {System.exit(0);}
		});

		window.setSize(250, 150);
		window.setVisible(true);


	}
}
