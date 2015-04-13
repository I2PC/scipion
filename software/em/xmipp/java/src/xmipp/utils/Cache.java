
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.utils;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author Juanjo Vega
 */
public class Cache<K, T> extends LinkedHashMap<K, T> {

    public final static int MEMORY_SIZE = 134217728;// 128 MB = 128 * 1024 * 1024 Bytes.
    public final static int MAXPXSIZE = 4;  // 4 Bytes => 4 * 8 = 32 bits as maximum pixel size.
    private int limit;  // Max number of elements in cache.

    public Cache() {
        this(1);
    }

    public Cache(int limit) {
        super(limit);

        this.limit = limit;
    }

    @SuppressWarnings("element-type-mismatch")
    public void resize(int limit) {
        // If new limit is lower than previous, remaining items are removed.
        DEBUG.printMessage(" *** Cache: resizing to: [" + limit + "] elements");
        //DEBUG.printMessage(" *** this.limit: [" + this.limit + "] / limit: [" + limit + "]");

        if (this.limit > limit && limit < size()) {
            Iterator iter = keySet().iterator();
            int items2remove = size() - limit;
            ArrayList toRemove = new ArrayList(items2remove);

            // To avoid synchronizing problems, store items in a temporary array...
            while (items2remove-- > 0) {
                toRemove.add(iter.next());
            }

            // ...and remove them later.
            for (Object o : toRemove) {
                remove(o);
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
        }

        return object;
    }

    /** When limit is reached, oldest item is removed */
    @Override
    protected boolean removeEldestEntry(Map.Entry eldest) {
        return size() > limit;
    }
}
