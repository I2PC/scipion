/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.windows.menubar;

import java.awt.Container;
import java.awt.Menu;
import java.awt.MenuBar;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.util.Vector;

/**
 *
 * @author Juanjo Vega
 */
public class DynamicMenuBar extends MenuBar {

    private final static int CHARWIDTH = 10;
    private final static int EXTENSIONWIDTH = CHARWIDTH * 4;
    private Container parent;
    private Vector<Menu> menus = new Vector<Menu>();
    // Extension menu to display menus that cannot be displayed in the menubar
    private Menu extension = new Menu(">>");

    public DynamicMenuBar(Container parent) {
        super();

        this.parent = parent;

        parent.addComponentListener(new ResizeListener());
    }

    public void addMenu(Menu m) {
        menus.add(m);

        updateGUI();
    }

    public void updateGUI() {
        DynamicMenuBar menu = DynamicMenuBar.this;

        // Clears menus.
        menu.removeAll();
        extension.removeAll();

        // Adds as much items as possible to the normal menu.
        int index = 0;
        while (index < menus.size()
                && menu.getWidth() + EXTENSIONWIDTH < parent.getWidth()) {
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
    }

    public void removeAll() {
        while (getMenuCount() > 0) {
            this.remove(0);
        }
    }

    public int getWidth() {
        int width = 0;

        for (int i = 0; i < getMenuCount(); i++) {
            width += menus.get(i).getLabel().length() * CHARWIDTH;
        }

        return width;
    }

    // Every time the menu is resized, check what menus can be displayed at their preferred size and move the rest to
    class ResizeListener extends ComponentAdapter {

        @Override
        public void componentResized(ComponentEvent evt) {
            updateGUI();
        }
    }
/*
    public static void main(String args[]) {
        String[] labels = {"File", "Edit", "Format", "View", "Tools", "Configure", "Admin", "Compile", "Options", "Help"};
        JFrame f = new JFrame("Test Layout");

        DynamicMenuBar mbar = new DynamicMenuBar(f);

        // Creates menus and items.
        for (int i = 0; i < labels.length; i++) {
            Menu menu = new Menu(i + "." + labels[i]);
            menu.add(new MenuItem("submenu_" + i));

            mbar.addMenu(menu);
        }

        f.setMenuBar(mbar);
        f.setSize(800, 200);
        f.addWindowListener(new WindowAdapter() {

            @Override
            public void windowClosing(WindowEvent e) {
                System.exit(0);
            }
        });
        f.setVisible(true);
    }*/
}
