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

// some recommend that plugin classes must be in the default package...
// here we use it has its own submenu
package xmipptomo;


import ij.*;
import ij.gui.*;
import ij.io.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.Calendar;
import ij.plugin.frame.*;

/**
 * @author jcuenca
 * extends PlugInFrame ImageJ plugin base class
 * implements ActionListener (buttons), MouseListener (mouse move & clicks)
 * 
 * Main class, responsible for workflow control and plugin initialization
 */

// underscore in the name is required for automatic installation in the plugins menu...
// Requirements: http://u759.curie.u-psud.fr/compteur/download.php?Fichier=software/update/20090928/U759_InputOutput.jar
public class Xmipp_Tomo extends PlugInFrame implements ActionListener,MouseListener {

	/**
	 * Required because of inherited "serializationability"...
	 */
	private static final long serialVersionUID = -4063711977454855701L;
	
	/** Protocol for adding new buttons to the workflow/UI:
	 *  - set the label in Commands enum as static String,
	 *  - run addButton in Xmipp_Tomo constructor, 
	 *  - update actionPerformed method
	 */
	
	// labels of buttons
	private static enum ButtonLabels {
		CMD_INFO("Info"),CMD_PREPROC("Preprocessing"), CMD_ALIGN("Align");
		private final String label;
		ButtonLabels(String s) { this.label=s;}
		public String label(){return label;}
	};
	
	private static enum CmdExitValues {
		OK(0),ERROR(1);
		private final int value;
		CmdExitValues(int err) { value=err;}
		public int value(){return value;}
	};
	
	// enable debuging tests - release versions should have this set to 0
	final static int TESTING = 1;
	
	// store all data following the MVC pattern
	private TomoData dataModel;

	// UI elements
	Panel panel; // main button panel
	// display tilt series
	TomoWindow tw=null;
	
	// implement UI concurrency with private class that extends Thread
	// ImportDataThread: load & display projections in parallel
	private class ImportDataThread extends Thread{
		private TomoData dataModel;
		private String dataPath;
		
		ImportDataThread(TomoData model,String path){
			dataModel=model;
			dataPath=path;
		}
		
		public void run() {
			try{
				dataModel.import_data(dataPath);
			}catch (IOException ex){
				debug("ImportDataThread.run - Error opening file");
			}catch (InterruptedException ex){
				debug("ImportDataThread.run - Interrupted exception");
			}catch (Exception ex){
				debug("ImportDataThread.run - unexpected exception", ex);
			}
		}
	}

	public Xmipp_Tomo() {
		super("Xmipp_Tomo");
		
		dataModel=new TomoData();
		
		setLayout(new FlowLayout());
		panel = new Panel();
		panel.setLayout(new GridLayout(4, 4, 5, 5));
		panel.add(new Label("Xmipp Tomo"));
		addButton(ButtonLabels.CMD_INFO.label());
		addButton(ButtonLabels.CMD_PREPROC.label());
		addButton(ButtonLabels.CMD_ALIGN.label());
		
		add(panel);
		
		pack();
		GUI.center(this);
		// show is deprecated...
		show();
	}
	
	/* (non-Javadoc)
	 * wait for user input (in button panel)
	 * @see ij.plugin.frame.PlugInFrame#run(java.lang.String)
	 */
	public void run(String arg){
	}

	
	/** Add button to main panel
	 * @param label
	 */
	void addButton(String label) {
		Button b = new Button(label);
		b.addActionListener(this);
		b.addKeyListener(IJ.getInstance());
		panel.add(b);
	}

	
	/** Run cmdline in a shell and show standard output & error via debug()
	 * @param cmdline
	 * @return the exit value of cmdline (@see Process.waitFor())
	 */
	int exec(String cmdline){
		// execution details may change with each OS...
		//String osName = System.getProperty("os.name" );
		
		int exitValue=CmdExitValues.OK.value();
		Process proc=null;
		
		Runtime rt = Runtime.getRuntime();
		
		try{
			proc = rt.exec(cmdline);
		}catch (IOException ex){
			return -1;
		}
		
		// Prepare buffered readers from inputstreams (stderr, stdout) ...
		InputStream stderr=proc.getErrorStream();
		InputStreamReader stderr_isr = new InputStreamReader(stderr);
        BufferedReader stderr_br = new BufferedReader(stderr_isr);

		InputStream stdout=proc.getInputStream();
		InputStreamReader stdout_isr = new InputStreamReader(stdout);
        BufferedReader stdout_br = new BufferedReader(stdout_isr);
		
		String line=null;
		
		try{
			// print stdout
			while ( (line = stdout_br.readLine()) != null)
	            debug(line);    
	        
			// print stderr
			while ( (line = stderr_br.readLine()) != null)
	            debug(line);    
			
        } catch (IOException ex){
            ex.printStackTrace();  
        }
        
        try{
        	exitValue = proc.waitFor();
        }catch (java.lang.InterruptedException ex){
        	exitValue=CmdExitValues.ERROR.value();;
        }
        return exitValue;
	} // exec end
	
	
	/** Show a file browser and... 
	 * @return the path of the file chosen by the user
	 */
	public String browseFile(){
		OpenDialog od = new OpenDialog("Import file",null);
        String directory = od.getDirectory();
		String fileName = od.getFileName();
		String path= directory + fileName;
		return path;
	}
	
	
	/* handle button/keyboard pressing, from both this plugin and the windows it opens
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	public void actionPerformed(ActionEvent e) {
	
		String label = e.getActionCommand();

		// select proper action method based on button's label
		if (label==null)
			return;
		else if (label.equals(ButtonLabels.CMD_INFO.label())){
			this.infoAction();
		}else if (label.equals(ButtonLabels.CMD_PREPROC.label())){
			exec("date");
		}else if (label.equals(ButtonLabels.CMD_ALIGN.label())){
			String directory=IJ.getDirectory("");
			exec("ls "+directory);
		}
	} // actionPerformed end
	
	
	/********************* Mouse event handlers ****************************/
	
	 public void mousePressed(MouseEvent e) {
	        // debug("Mouse pressed; # of clicks: "  + e.getClickCount());
	 }

	
	public void mouseReleased(MouseEvent e) {
	}

	    public void mouseEntered(MouseEvent e) {
    }

	    public void mouseExited(MouseEvent e) {
	}

    public void mouseClicked(MouseEvent e) {
	}

    /********************* Trace methods for debugging ****************************/
    
    
	/** print s with a timestamp, using IJ Results window (@see IJ.write)
	 * @param s
	 */
	public static void debug(String s){
		Calendar calendar = Calendar.getInstance();
		java.util.Date now = calendar.getTime();
		IJ.write("" + now.getTime() + " > " + s);
	}

	/** same as Debug, plus ex stack trace
	 * @param s
	 * @param ex Exception
	 */
	public static void debug(String s, Exception ex){
		debug(s);
		ex.printStackTrace();
	}
	
	public TomoData getDataModel() {
		return dataModel;
	}

	public void setDataModel(TomoData dataModel) {
		this.dataModel = dataModel;
	}

	
	private void createTomoWindow(){
		if(getDataModel().getImage()==null)
			debug("No image available to display");
		else{	
			if(tw==null){
				tw= new TomoWindow();
				tw.create(getDataModel());
				WindowManager.addWindow(tw);
			}
		}
	}
		
	
	/* Procedures handling application workflows (in response to events) */
	

	/**
	 * All actions corresponding to the workflow "Info"
	 */
	private void infoAction(){
		String path = browseFile();
		
		try{
			// import data in one thread and hold this thread until the first projection is loaded
			(new Thread(new ImportDataThread(getDataModel(),path))).start();
			getDataModel().waitForFirstImage();
		}catch (InterruptedException ex){
			debug("Xmipp_Tomo - Interrupted exception");
		}

		createTomoWindow();
		if(tw != null)
			tw.display();
	}
	
	// ij.IJ class offers a lot of useful methods...
	private void ijExamples(){
	
		IJ.showStatus("Xmipp_Tomo started.");
		IJ.showProgress(0.0);
		String name = IJ.getString("Please enter your name: ","I.J. User");
		IJ.showProgress(0.5);
		IJ.write("Starting sample plugin Red And Blue ... ");
		IJ.runPlugIn("Red_And_Blue","");
		IJ.showProgress(1.0);
		IJ.showMessage("Finished.",name+", thank you for running this plugin");
		
	}
	
}

