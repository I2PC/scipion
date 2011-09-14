
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author Juanjo Vega
 */
public class Cache<K, T> extends LinkedHashMap<K, T> {

    /** Max number of elements in cache */
    private int limit = 50;

    @SuppressWarnings("element-type-mismatch")
    public void resize(int limit) {
        // If new limit is lower than previous, remaining items are removed.
        Logger.debug(" *** resizing to: [" + limit + "]");
        Logger.debug(" *** this.limit: [" + this.limit + "] / limit: [" + limit + "]");

        if (this.limit > limit && limit < size()) {
            int i = 0;
            Iterator iter = this.keySet().iterator();
            int items2remove = size() - limit;

            Object trash[] = new Object[items2remove];

            while (iter.hasNext() & i < items2remove) {
                trash[i++] = iter.next();
            }

            for (int j = 0; j < trash.length; j++) {
                Object object = trash[j];
                remove(object);
            }
        }

        this.limit = limit;
    }

    @Override
    public T get(Object key) {
        T object = super.get((K) key);

        if (object != null) {   // Object has been referenced, so place it at higher position.
            remove((K) key);
            put((K) key, (T) object);
            Logger.debug("cache hit: " + key.toString());
        }

        return object;
    }

    /** When limit is reached, oldest item is removed */
    @Override
    protected boolean removeEldestEntry(Map.Entry eldest) {
        return size() > limit;
    }
}
