/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.coss;

/**
 *
 * @author Juanjo Vega
 */
public class CossEntry {

    protected boolean enabled;
    protected String img1, img2, img3;
    protected String info;

    public CossEntry(boolean enabled, String img1, String img2, String img3, String info) {
        this.enabled = enabled;
        this.img1 = img1;
        this.img2 = img2;
        this.img3 = img3;
        this.info = info;
    }
}
