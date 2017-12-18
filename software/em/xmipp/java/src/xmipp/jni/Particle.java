package xmipp.jni;

import java.awt.Point;



public class Particle implements Comparable<Particle> {
	protected int x;
	protected int y;
	protected double cost;
	
	public Particle(int x, int y)
	{
		this.x = x;
		this.y = y;
	}
        
        public Particle(int x, int y, double cost)
	{
		this.x = x;
		this.y = y;
		this.cost = cost;
	}
	
	
	public int getX() {
		return x;
	}

	public void setX(int x) {
		this.x = x;
	}

	public int getY() {
		return y;
	}

	public void setY(int y) {
		this.y = y;
	}
        
        public double getCost()
	{
		return cost;
	}
        
        public void setCost(double cost)
	{
		this.cost = cost;
		
	}

	
		
	public boolean contains(int x2, int y2, int size )
	{
		int radius = size/2;
			if(x2 < x - radius || x2 > x + radius)
				return false;
			if(y2 < y - radius || y2 > y + radius)
				return false;
			return true;
	}

	public void setPosition(int x, int y) {
		this.x = x;
		this.y = y;
	}
	
	public Point getPosition() {
		return new Point(x, y);
		
	}
	
	public String toString()
	{
		return String.format("x = %s; y = %s cost = %.2f", x, y, cost); 
	}
	

	@Override
	public int compareTo(Particle p) {
		if(p.x > x)
			return 1;
		if(p.x == x)
		{
			if(p.y > y)
				return 1;
			if(p.y == y)
				return 0;
			return -1;
		}
		return -1;
	}
}
