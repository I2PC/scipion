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
 * - Why?
 * It acts like a sweet bridge between Swing components events and our Controller.
 * Just create a TomoAction and pass it to the component's constructor
 * @extends AbstractAction - so you can link a TomoAction to a Swing component (say, a JButton). When the 
 *   component fires an event, TomoAction.actionPerformed is called
 * Alternatives to "Action paradigm": trying to identify a component in a "if list" by its label
 */

public class Action extends AbstractAction {
	// Following MVC paradigm, a mechanism is needed to establish which Controller's method to call when
	// some specific event (like clicking on a button) happens, and hence actionPerformed is called.
	// I opted for Java reflection, which allows for calling methods using a String.
	// The alternative would be a long "if list", again
	private Method methodToCall;
	private Object controller;
	// right now it's not necessary because controller's method have no
	// parameters
	protected Object[] params=null;

	// no need for the controller to be a TomoControler...
	// TODO: but maybe Object is too generic... use AbstractController?
	/**
	 * @param c - controller class which has the method that will be
	 * called in actionPerformed
	 */
	Action(Object c, Command cmd) {
		this(c,cmd,null);	
	}

	Action(Object c, Command cmd,Object parameter) {
		super(cmd.getLabel(), cmd.getIcon());
		setEnabled(cmd.isEnabled());
		controller = c;
		setMethodToCall(cmd, parameter);
		if(parameter != null){
			params=new Object[1];
			params[0]=parameter;
		}
	}

	private void setMethodToCall(Command cmd, Object parameter){
		try {
			if(parameter==null)
				methodToCall = controller.getClass().getMethod(cmd.getMethod());
			else
				methodToCall = controller.getClass().getMethod(cmd.getMethod(),parameter.getClass());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/*
	 * Calls the controller's method corresponding to this action
	 * @see
	 * java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	public void actionPerformed(ActionEvent ae) {
		try {
			if(params==null)
				methodToCall.invoke(controller);
			else
				methodToCall.invoke(controller, params);
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		} catch (NullPointerException ex){
			Logger.debug("TomoAction.actionPerformed " + ex.toString());
		}
	}

}
