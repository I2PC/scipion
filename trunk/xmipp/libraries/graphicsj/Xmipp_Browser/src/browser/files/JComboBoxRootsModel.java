/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.files;

import java.io.File;
import javax.swing.DefaultComboBoxModel;

/**
 *
 * @author Juanjo Vega
 */
public class JComboBoxRootsModel extends DefaultComboBoxModel {

    protected File roots[];

    public JComboBoxRootsModel() {
        update();
    }

    public void update() {
        roots = File.listRoots();
    }

    @Override
    public Object getElementAt(int index) {
        return roots[index];
    }

    @Override
    public int getSize() {
        return roots.length;
    }
}
