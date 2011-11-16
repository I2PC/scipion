package sphere;

import java.util.ArrayList;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class SphereCell {

    public ArrayList<String> fileNames;    // Image file name.
    public double value = 0.0;
    public boolean locked = false;

    public SphereCell() {
        fileNames = new ArrayList<String>();
    }
}
