/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tests.awtmenubar;

import java.awt.Container;
import java.awt.Menu;
import java.awt.MenuBar;
import java.awt.MenuItem;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.util.Vector;

/**
 *
 * @author Juanjo Vega
 */
public class DynamicJMenuBar extends MenuBar {

    private Container parent;
    private Vector<Menu> menus = new Vector<Menu>();
    // Extension menu to display menus that cannot be displayed in the menubar
    private Menu extension = new Menu(">>");

    public DynamicJMenuBar(Container parent, String items[]) {
        super();

        this.parent = parent;

        // Creates menus and items.
        for (int i = 0; i < items.length; i++) {
            Menu menu = new Menu(i + "." + items[i]);
            menu.add(new MenuItem("submenu_" + i));

            menus.add(menu);
        }

//        setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));
//        addComponentListener(new ResizeListener());
        parent.addComponentListener(new ResizeListener());
    }

    public void removeAll() {
        while (getMenuCount() > 0) {
            this.remove(0);
        }
    }

    // Every time the menu is resized, check what menus can be displayed at their preferred size and move the rest to
    class ResizeListener extends ComponentAdapter {

        @Override
        public void componentResized(ComponentEvent evt) {
            DynamicJMenuBar menu = DynamicJMenuBar.this;

            // Clears menus.
            menu.removeAll();
            extension.removeAll();

            // Adds as much items as possible to the normal menu.
            int index = 0;
            /*
            while (index < menus.size()
            && menu.getPreferredSize().width + extension.getPreferredSize().width
            < menu.getWidth()) {
             */
            while (index < menus.size()
                    && menu.getMenuCount() * 60 + 10//extension.getPreferredSize().width
                    < parent.getWidth()) {

                menu.add(menus.get(index++));
                System.out.println(" MENU ADDED: " + menus.get(index - 1).getLabel());
            }

            if (index != menus.size()) {
                if (index > 0) {
                    menu.remove(--index);
                }
                menu.add(extension);
            }

            while (index < menus.size()) {
                extension.add(menus.get(index++));
                System.out.println(" EXTENSION [!] ITEM ADDED: " + menus.get(index - 1).getLabel());
            }

            System.out.println(" ----------------------------------------------- ");
            //menu.updateUI();
        }
    }
}
