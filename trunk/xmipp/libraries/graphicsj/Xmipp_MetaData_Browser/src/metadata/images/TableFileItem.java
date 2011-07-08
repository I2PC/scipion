/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package metadata.images;

/**
 *
 * @author Juanjo Vega
 */
public class TableFileItem {

    protected String path;
    protected String originalValue;

    public TableFileItem(String path, String originalValue) {
        this.path = path;
        this.originalValue = originalValue;
    }

    public String getPath() {
        return path;
    }

    public String getOriginalValue() {
        return originalValue;
    }
}
