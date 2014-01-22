/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
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

package xmipp.utils;

import java.awt.Dialog;
import java.awt.event.ActionEvent;

import javax.swing.JButton;
import javax.swing.JFrame;

/** Special type of message with Yes, No, Cancel options */
public class XmippQuestionDialog extends XmippMessageDialog {
	private static final long serialVersionUID = 1L;
	
	protected JButton btnQuestionCancel; //Button cancel will be 'No' button
	protected String btnQuestionCancelText = "Cancel";
	protected boolean btnQuestionCancelDisplay = true;
	protected boolean canceled = false;
	private static final String TITLE = "QUESTION";
	private static final String ICON = "question.gif";
        
        public static int YES_OPTION = 0;
        public static int NO_OPTION = 1;
        public static int CANCEL_OPTION = 2;
        protected int option;
	
	public XmippQuestionDialog(JFrame parent, String message) {
		super(parent, TITLE);
		init(message, ICON);	
	}
	
	public XmippQuestionDialog(Dialog parent, String message) {
		super(parent, TITLE);
		init(message, ICON);	
	}
	
	public XmippQuestionDialog(JFrame parent, String message, boolean showCancel) {
		super(parent, TITLE);
		btnQuestionCancelDisplay = showCancel;
		init(message, ICON);	
	}
	
	public XmippQuestionDialog(Dialog parent, String message, boolean showCancel) {
		super(parent, TITLE);
		btnQuestionCancelDisplay = showCancel;
		init(message, ICON);		
	}
	
	@Override
	public void handleActionPerformed(ActionEvent evt) {
            JButton button = (JButton)evt.getSource();
		if (btnQuestionCancel == button){
                    canceled = true;
                    close(false);
		}
                
	}
	
	@Override
	protected void createButtons() {
		btnOkText = "Yes";
		btnCancelText = "No";
		super.createButtons();
		if (btnQuestionCancelDisplay) {
			btnQuestionCancel = XmippWindowUtil.getTextButton(btnQuestionCancelText, this);
			panelBtn.add(btnQuestionCancel);
                        btnQuestionCancel.addActionListener(this);
		}
	}
	
	public boolean isCanceled(){
		return canceled;
	}
        
        public int getOption(){
            if (result)
                return YES_OPTION;
            if (isCanceled())
                return CANCEL_OPTION;
            return NO_OPTION;
        }
}//class XmippQuestionDialog