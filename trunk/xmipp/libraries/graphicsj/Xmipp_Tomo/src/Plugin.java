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
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;


/**
 * @author jcuenca
 * Interface to ImageJ plugins - collects parameters from IJ plugin dialogs to run the plugin with them later
 * Each subclass knows the specific parameters of a plugin
 * 
 * Protocol for adding a plugin
 * - Create a subclass of Plugin (see GaussianPlugin for example)
 * - Add the command definition in Command
 * - add the action method in TomoController
 * - Add a button (see TomoWindow)
 */
public abstract class Plugin {
	// TODO: background plane substraction (@see Xmipp filters.h)
	// TODO: Roman filter
	
	/** Analyze the dialog gd and save its parameters
	 *  Note: when using the gd.getNextXX() methods, if the user simply presses OK, the method returns 0
	 *  (instead of the default value)
	 * @param gd
	 */
	public abstract void collectParameters(GenericDialog gd);
	
/**
 * The name of the command is the text exactly as it is displayed in the IJ Menu (without including
 * menus and submenus, just the command label)
 * @return
 */
	public abstract String getCommand();
	/**
	 * Note: the names of the parameters (required to set macro options) are the labels of the text fields associated to them.
	 *       Syntax: param1=value1 param2=value2 (separator is blank space)
	 *       @see GaussianPlugin, Macro.getValue
	 * @return a string with the plugin options in a format that IJ.run understands (Macro format)
	 */
	public abstract String getOptions();
	

	
	public final void run(ImagePlus image){
		WindowManager.setTempCurrentImage(image);
		WindowManager.setWindow(null);
		IJ.run(getCommand(),getOptions());
	}
	

}
