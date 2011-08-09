package sphere;




import java.util.Vector;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class SphereCell {

    public Vector<String> fileNames;    // Image file name.
    public double value = 0.0;
    public boolean locked = false;

    public SphereCell() {
        fileNames = new Vector<String>();
    }
}
