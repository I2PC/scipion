/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.tableitems;

/**
 *
 * @author Juanjo Vega
 */
public interface iMDComparable<T> {//extends Comparable<T> {

    public int compareToByLabel(T object, int label);
}
