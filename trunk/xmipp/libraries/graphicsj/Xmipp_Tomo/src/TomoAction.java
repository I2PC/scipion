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

import java.awt.event.ActionEvent;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import javax.swing.AbstractAction;

/**
 * Handles user interactions connecting commands with the window controller, with the help of Java reflection
 * 
 */

public class TomoAction extends AbstractAction {
	private Method methodToCall;
	private TomoController controller;
	// right now it's not necessary because controller's method have no
	// parameters
	protected Object[] params;

	TomoAction(TomoController c, Command cmd) throws Exception{
		super(cmd.getLabel(), cmd.getIcon());
		setEnabled(cmd.isEnabled());
		controller = c;
	    try {
	      methodToCall = controller.getClass().getMethod(cmd.getMethod());
	    } catch (Exception e) {
	      e.printStackTrace();
	    }
	}

	/*
	 * Calls a controller's method according to the command cmd
	 * 
	 * @see
	 * java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	public void actionPerformed(ActionEvent ae) {
		try {
			methodToCall.invoke(controller);
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		} catch (NullPointerException ex){
			Xmipp_Tomo.debug("TomoAction.actionPerformed " + ex.toString());
		}
	}

}
