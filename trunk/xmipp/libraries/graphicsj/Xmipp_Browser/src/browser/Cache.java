
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser;

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
//        System.out.println(" *** resizing to: [" + limit + "]");
//        System.out.println(" *** this.limit: [" + this.limit + "] / limit: [" + limit + "]");
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
        }

        return object;
    }

    /** When limit is reached, oldest item is removed */
    @Override
    protected boolean removeEldestEntry(Map.Entry eldest) {
        return size() > limit;
    }
    /*
    static void traverse(Collection coll) {
    Iterator iter = coll.iterator();
    while (iter.hasNext()) {
    System.out.print(iter.next() + " ");
    }
    System.out.println();
    }

    public static void main(String args[]) {
    Cache cache = new Cache(5);

    for (int i = 0; i < 5; i++) {
    cache.put(String.valueOf(i), i);
    }
    traverse(cache.entrySet());

    System.out.println("elem[1]: " + cache.get(String.valueOf(1)));
    traverse(cache.entrySet());
    //cache.resize(3);
    cache.put(5, 5);
    traverse(cache.entrySet());
    cache.put(6, 6);
    traverse(cache.entrySet());
    cache.put(7, 7);

    traverse(cache.entrySet());

    System.out.println("elem[4]: " + cache.get(String.valueOf(4)));
    traverse(cache.entrySet());
    cache.resize(3);
    traverse(cache.entrySet());
    /*
    cache.put(11, 11);
    cache.traverse(cache.entrySet());
    cache.get(3);
    cache.traverse(cache.entrySet());
    cache.put(12, 12);

    cache.traverse(cache.entrySet());
    cache.resize(4);

    cache.traverse(cache.entrySet());
    }*/
}
