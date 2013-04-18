package xmipp.viewer.ctf;

import ij.ImagePlus;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

import javax.swing.JPanel;

public class CTFAnalyzerImagePane extends JPanel implements MouseListener, MouseWheelListener, MouseMotionListener
{

	private int x0;
	private int y0;
	private int x1;
	private int y1;
	private double profileangle;
	private ImagePlus imp;
	private int radius;
	private CTFAnalyzerJFrame frame;


	public CTFAnalyzerImagePane(ImagePlus imp, CTFAnalyzerJFrame frame)
	{
		this.imp = imp;
		this.frame = frame;
		setPreferredSize(new Dimension(imp.getWidth(), imp.getHeight()));

		radius = imp.getWidth() / 2;
		x0 = radius;
		y0 = radius;
		profileangle = 0;
		x1 = x0 + (int) (radius * Math.cos(profileangle));
		y1 = y0 + (int) (radius * Math.sin(profileangle));
		setToolTipText("Use mouse to move profile");
		addMouseListener(this);
		addMouseWheelListener(this);
		addMouseMotionListener(this);
	}
	
	
	public double getProfileangle()
	{
		return profileangle;
	}


	public int getX0()
	{
		return x0;
	}


	public int getY0()
	{
		return y0;
	}


	public int getX1()
	{
		return x1;
	}


	public int getY1()
	{
		return y1;
	}


	public void paintComponent(Graphics g)
	{
		super.paintComponent(g);
		
		g.drawImage(imp.getImage(), 0, 0, this);
		g.setColor(Color.green);
		g.drawLine(x0, y0, x1, y1);
		
	}
	
	private void updateProfile(MouseEvent e)
	{
		int y = radius - e.getY();
		int x = e.getX() - radius;
		double angle = -Math.atan2(y, x);
		updateProfile(angle);
	}

	private void updateProfile(double angle)
	{
		profileangle = angle;
		x1 = x0 + (int) (radius * Math.cos(profileangle));
		y1 = y0 + (int) (radius * Math.sin(profileangle));
		repaint();
		frame.fillGraphics();
	}


	@Override
	public void mouseWheelMoved(MouseWheelEvent e)
	{
		int rotation = e.getWheelRotation();//1 or -1
		updateProfile(profileangle + rotation * Math.toRadians(2));
		
	}


	@Override
	public void mouseClicked(MouseEvent e)
	{
		updateProfile(e);
		
	}
	
	


	@Override
	public void mouseEntered(MouseEvent arg0)
	{
		// TODO Auto-generated method stub
		
	}


	@Override
	public void mouseExited(MouseEvent arg0)
	{
		// TODO Auto-generated method stub
		
	}


	@Override
	public void mousePressed(MouseEvent arg0)
	{
		// TODO Auto-generated method stub
		
	}


	@Override
	public void mouseReleased(MouseEvent e)
	{
		updateProfile(e);
		
	}


	@Override
	public void mouseDragged(MouseEvent e)
	{
		updateProfile(e);
		
	}


	@Override
	public void mouseMoved(MouseEvent arg0)
	{
		// TODO Auto-generated method stub
		
	}

}
