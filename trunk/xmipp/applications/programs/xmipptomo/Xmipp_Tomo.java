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
import java.util.Date;




import ij.plugin.frame.*;


// underscore in the name is required for automatic installation in the plugins menu...
// Requirements: http://u759.curie.u-psud.fr/compteur/download.php?Fichier=software/update/20090928/U759_InputOutput.jar

public class Xmipp_Tomo extends PlugInFrame implements ActionListener,MouseListener {

	/**
	 * Required because of inherited serializationability...
	 */
	private static final long serialVersionUID = -4063711977454855701L;
	/** Adding buttons: set the label here as static String, addButton in constructor, 
	 *  update actionPerformed method
	 */
	private static enum Commands {
		CMD_INFO("Info"),CMD_PREPROC("Preprocessing"), CMD_ALIGN("Align");
		private final String label;
		Commands(String s) { this.label=s;}
		public String label(){return label;}
	};
	
	final static int PLUGIN_NOT_FOUND = -1;
	// enable debuging tests - release versions should have this set to 0
	final static int TESTING = 1;
	
	private TomoData dataModel;

	// UI elements
	Panel panel; // main button panel
	
	TomoWindow tw=null;
	
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
				debug("ImportDataThread@Xmipp_Tomo - Error opening file");
			}catch (InterruptedException ex){
				debug("ImportDataThread@Xmipp_Tomo - Interrupted exception");
			}catch (Exception ex){
				debug("ImportDataThread@Xmipp_Tomo - unexpected exception", ex);
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
		addButton(Commands.CMD_INFO.label());
		addButton(Commands.CMD_PREPROC.label());
		addButton(Commands.CMD_ALIGN.label());
		
		add(panel);
		
		pack();
		GUI.center(this);
		// show is deprecated...
		show();
	}
	
	// wait for user input (in button panel)
	public void run(String arg) {
		// ij.IJ class offers all these useful methods.
		/* IJ.showStatus("Xmipp_Tomo started.");
		IJ.showProgress(0.0);
		
		String name = IJ.getString("Please enter your name: ","I.J. User");
		IJ.showProgress(0.5);
		IJ.write("Starting sample plugin Red And Blue ... ");
		IJ.runPlugIn("Red_And_Blue","");
		IJ.showProgress(1.0);
		IJ.showMessage("Finished.",name+", thank you for running this plugin");
		*/
		
	}

	void addButton(String label) {
		Button b = new Button(label);
		b.addActionListener(this);
		b.addKeyListener(IJ.getInstance());
		panel.add(b);
	}
	
	int exec(String cmdline){
		// execution details may change with each OS...
		//String osName = System.getProperty("os.name" );
		int exitValue=0;
		Process proc=null;
		
		Runtime rt = Runtime.getRuntime();
		
		try{
			proc = rt.exec(cmdline);
		}catch (IOException ex){
			return -1;
		}
		// Prepare buffered readers from inputstreams...
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
        	exitValue=-1;
        }
        return exitValue;
	} // exec end
	
	public String browseFile(){
		OpenDialog od = new OpenDialog("Import file",null);
        String directory = od.getDirectory();
		String fileName = od.getFileName();
		String path= directory + fileName;
		// IJ.write(path);
		return path;
	}
	

	// handle button/keyboard pressing, from both this plugin and the windows it opens
	public void actionPerformed(ActionEvent e) {
	
		String label = e.getActionCommand();
		// debug(label);
		if (label==null)
			return;
		else if (label.equals(Commands.CMD_INFO.label())){
			this.infoAction();
			
			// IJ.showMessage("Running "+CMD_INFO);
			//import_data();		
			//getImage().show();
			/* getImage().getWindow().addMouseListener(this);
			Component [] cs=getImage().getWindow().getComponents();
			for(int i=0;i<cs.length; i++)
				cs[i].addMouseListener(this); */
		}else if (label.equals(Commands.CMD_PREPROC.label())){
			// IJ.showMessage("Running "+CMD_PREPROC);
			exec("date");
		}else if (label.equals(Commands.CMD_ALIGN.label())){
			String directory=IJ.getDirectory("");
			exec("ls "+directory);
		}
	} // actionPerformed end
	
	 public void mousePressed(MouseEvent e) {
	        // debug("Mouse pressed; # of clicks: "  + e.getClickCount());
	 }

	
	public void mouseReleased(MouseEvent e) {
	     // debug("Mouse released; # of clicks: " + e.getClickCount());
	}

	    public void mouseEntered(MouseEvent e) {
	       // debug("Mouse entered", e);
    }

	    public void mouseExited(MouseEvent e) {
	       // debug("Mouse exited", e);
	}

	    public void mouseClicked(MouseEvent e) {
	    	// String componentName=e.getComponent().getName();
	    	// if("scrollbar0".equals(componentName))
	    	//	debug("Mouse clicked (# of clicks: " + e.getClickCount() + componentName +")");
	    	//	getImage().getWindow().getGraphics().drawString("KK", 5, 5);
	}

	public static void debug(String s){
		Calendar calendar = Calendar.getInstance();
		java.util.Date now = calendar.getTime();
		IJ.write("" + now.getTime() + " > " + s);
	}

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
	
	/* All actions corresponding to the workflow "Info" */
	private void infoAction(){
		String path = browseFile();
		
		try{
			(new Thread(new ImportDataThread(getDataModel(),path))).start();
			//getDataModel().import_data(path);
			getDataModel().waitForFirstImage();
		}catch (InterruptedException ex){
			debug("Xmipp_Tomo - Interrupted exception");
		}

		createTomoWindow();
		if(tw != null)
			tw.display();
	}
	
}

