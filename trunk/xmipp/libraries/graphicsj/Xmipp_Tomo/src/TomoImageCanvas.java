/***************************************************************************
 * 
 * @author: Jesus Cuenca (jcuenca@cnb.csic.es)
 * 
 *          Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 * 
 *          This program is free software; you can redistribute it and/or modify
 *          it under the terms of the GNU General Public License as published by
 *          the Free Software Foundation; either version 2 of the License, or
 *          (at your option) any later version.
 * 
 *          This program is distributed in the hope that it will be useful, but
 *          WITHOUT ANY WARRANTY; without even the implied warranty of
 *          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *          General Public License for more details.
 * 
 *          You should have received a copy of the GNU General Public License
 *          along with this program; if not, write to the Free Software
 *          Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 *          USA
 * 
 *          All comments concerning this program package may be sent to the
 *          e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
/**
 * - Why?
 * XmippTomo needs some visual customizations which only worked with a custom Canvas
 */

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Line2D;
import ij.gui.ImageCanvas;

public class TomoImageCanvas extends ImageCanvas {
	private StackModel model;

	public TomoImageCanvas(StackModel model) {
		super(model.getCurrentImage());
		this.model = model;
	}

	public void paint(Graphics g) {
		super.paint(g);

		if (model != null) {
			if (model.isCurrentEnabled() == false) {
				drawDiscardedMark(g);
			}
		}
	}

	/**
	 * Overlay a cross (X) in the canvas
	 * @param g
	 */
	// TODO: keep the ROI displayed in the StackViewer when the user changes the slice (for instance, clicking "Play")
	private void drawDiscardedMark(Graphics g){
		Graphics2D g2 = (Graphics2D) g;
		g2.setColor(Color.RED);
		g2.setStroke(new BasicStroke(4));
		Line2D line1 = new Line2D.Float(5, 5, getWidth() - 5, getHeight() - 5);
		Line2D line2 = new Line2D.Float(getWidth() - 5, 5, 5, getHeight() - 5);
		g2.draw(line1);
		g2.draw(line2);
	}
}
