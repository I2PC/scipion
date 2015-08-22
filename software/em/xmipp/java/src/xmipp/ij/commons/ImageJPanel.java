/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.ij.commons;

import ij.ImagePlus;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import javax.swing.JPanel;

/**
 *
 * @author airen
 */
public class ImageJPanel extends JPanel{

    protected ImagePlus imp;
    protected int width, height;

    public ImageJPanel(ImagePlus image, int width, int height) {
       this.imp = image;
       this.width = width;
       this.height = height;
    }
    
    @Override
    public Dimension getPreferredSize()
    {
        return new Dimension(width, height);
        
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        
        g.drawImage(imp.getImage().getScaledInstance(width, height, Image.SCALE_SMOOTH), 0, 0, null); // see javadoc for more info on the parameters            
    }
    
    public void setMask(ImagePlus mask, int width, int height)
    {
        this.imp = mask;
        this.width = width;
        this.height = height;
        setSize(new Dimension(width, height));
        repaint();
        
    }

    

}