
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author Juanjo Vega
 */
public class Cache<K, T> extends LinkedHashMap<K, T> {

    /** Max number of elements in cache */
    private int limit;

    public Cache(int limit) {//, boolean show) {
        super();

        this.limit = limit;
    }

    public void setLimit(int limit) {
//        System.out.println("L=" + limit);
        this.limit = limit;
    }

    /** When limit is reached, oldest item is removed */
    @Override
    protected boolean removeEldestEntry(Map.Entry eldest) {
        return size() > this.limit;
    }
}
