/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tests;

import java.awt.FlowLayout;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.util.Vector;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;

/**
 *
 * @author Juanjo Vega
 */
public class DynamicJMenuBar extends JMenuBar {

    private Vector<JMenu> menus = new Vector<JMenu>();
    // Extension menu to display menus that cannot be displayed in the menubar
    private JMenu extension = new JMenu(">>");

    public DynamicJMenuBar(String items[]) {
        super();

        // Creates menus and items.
        for (int i = 0; i < items.length; i++) {
            JMenu menu = new JMenu(i + "." + items[i]);
            menu.add(new JMenuItem("submenu_" + i));

            menus.add(menu);
        }

        setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));
        addComponentListener(new ResizeListener());
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
            while (index < menus.size()
                    && menu.getPreferredSize().width + extension.getPreferredSize().width
                    < menu.getWidth()) {

                menu.add(menus.get(index++));
            }

            if (index != menus.size()) {
                if (index > 0) {
                    menu.remove(--index);
                }
                menu.add(extension);
            }

            while (index < menus.size()) {
                extension.add(menus.get(index++));
            }

            menu.updateUI();
        }
    }
}
