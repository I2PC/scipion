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

import java.io.File;
import javax.swing.JFileChooser;

import xmipp.jni.Filename;


/** Xmipp customization of the default XmippFileChooser */
@SuppressWarnings("serial")
public class XmippFileChooser extends JFileChooser {
	
	/** Constructor */	
	public XmippFileChooser(){
		this(new File(Filename.getCurrentDir()));
	}
	
	public XmippFileChooser(File f){
		super(f);
		init();
	}

    public XmippFileChooser(String path) {
        this(new File(path));
    }
	
	/** Perform customization */
	private void init(){		
		//Add custom icons for file types.
		setFileView(new XmippFileView());
		
		//Add the preview pane.
		setAccessory(new XmippFilePreview(this));
	}
	
	/** Shortcut function to return selected path */
	public String getSelectedPath(){
		return getSelectedFile().getPath();
	}

}
