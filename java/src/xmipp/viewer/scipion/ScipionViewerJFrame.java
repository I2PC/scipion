/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.awt.Dimension;
import java.awt.Image;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import xmipp.utils.ImageRenderer;
import xmipp.utils.XmippWindowUtil;


/**
 *
 * @author airen
 */
class ScipionViewerJFrame extends JFrame{
    
    protected ScipionMetaData md;
    protected JTable table;
    private ScipionViewerTableModel model;

    public ScipionViewerJFrame(ScipionMetaData md) {
        this.md = md;
        initComponents();
    }

    protected void initComponents() {
        setTitle(md.getDBFile());
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
       
        model = new ScipionViewerTableModel(md);
        table = new JTable(model);
        table.setDefaultRenderer(ImageIcon.class, new ImageRenderer());
        table.setRowHeight(md.getCellHeight());
        for(int i = 0; i < table.getColumnCount(); i ++)
            table.getColumnModel().getColumn(i).setPreferredWidth(md.getCellWidth());
        table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        JScrollPane sp = new JScrollPane(table);
        sp.setViewportView(table);
        sp.setPreferredSize(new Dimension(800, 600));
        add(sp);
        pack();
        XmippWindowUtil.setLocation(0.5, 0.2f, this);
        setVisible(true);
    }
    
}
