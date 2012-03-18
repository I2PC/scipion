/***************************************************************************
 * Authors:     Juanjo Vega
 * 				J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

package xmipp.viewer;

import java.awt.Image;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author Juanjo Vega
 */
public class RotSpectraItem extends ImageItem{

    JFreeChart chart;
    double vector[];

    public RotSpectraItem(String k, String l, double vector[]) {
        super(k, l);
        this.vector = vector;
        chart = createChart(l, vector);
    }

    /** Create chart based on vector data */
    static JFreeChart createChart(String label, double vector[]) {
        XYSeries series = new XYSeries(label);

        for (int i = 0; i < vector.length; i++) 
            series.add((double) i, vector[i]);

        XYSeriesCollection dataset = new XYSeriesCollection(series);

        JFreeChart chart = ChartFactory.createXYLineChart(
                "", "", "",
                dataset, PlotOrientation.VERTICAL,
                true, true, false);

        chart.removeLegend();

        return chart;
    }//function createChart

    @Override
    public Image getImage() {
        return chart.createBufferedImage(cellDim.width, cellDim.height);
    }

    public JFrame getChart() {
        ChartPanel panel = new ChartPanel(chart);
        JFrame frame = new JFrame();
        frame.setTitle(label);
        frame.getContentPane().add(panel);
        frame.pack();
        return frame;
    }
}//class RotSpectraItem
