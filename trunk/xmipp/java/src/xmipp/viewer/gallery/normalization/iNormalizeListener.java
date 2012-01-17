/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.gallery.normalization;

/**
 *
 * @author Juanjo Vega
 */
public interface iNormalizeListener {

    public void normalize(double min, double max);

    public void setNormalizedAuto();

    public void disableNormalization();
}
