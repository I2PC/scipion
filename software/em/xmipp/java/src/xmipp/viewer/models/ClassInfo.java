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

package xmipp.viewer.models;

import java.awt.Color;
import javax.swing.Icon;
import xmipp.utils.ColorIcon;

/** Structure to store class info */
public class ClassInfo {
	public String comment;
	public ColorIcon icon;
	public int index; // index of the class
	public int numberOfClasses; // Classes assigned to superclass
	public long numberOfImages; // total images assigned to superclass

	/** Constructor */
	public ClassInfo(String name, Color c) {
		comment = name;
		icon = new ColorIcon(c, 16, 16, 3, true, true);
	}
        
        /** Constructor */
	public ClassInfo(String name, Color c, int numberOfImages) {
		comment = name;
		icon = new ColorIcon(c, 16, 16, 3, true, true);
                this.numberOfImages = numberOfImages;
	}

	public String getComment() {
		return comment;
	}

	public void setComment(String value) {
		comment = value;
	}

	public Color getColor() {
		return icon.getColor();
	}

	public void setColor(Color value) {
		icon.setColor(value);
	}

	public Icon getIcon() {
		return icon;
	}
}// class ClassInfo