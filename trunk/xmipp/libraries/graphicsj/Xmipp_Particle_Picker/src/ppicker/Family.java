package ppicker;
import java.awt.Color;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;


public class Family {
	
	private String name;
	private Color color;
	private int radius;
	private static Family dfamily = new Family("Any", Color.green);
	private static Color[] colors = new Color[]{Color.BLUE, Color.CYAN, 
										Color.GREEN,
										Color.MAGENTA, Color.ORANGE, 
										Color.PINK, Color.RED, Color.YELLOW};
	private List<Particle> particles;
	
	public Family(String name, Color color, int radius)
	{
		this.name = name;
		this.color = color;
		this.radius = radius;
		particles = new ArrayList<Particle>();
	}
	
	public Family(String name, Color color)
	{
		this.name = name;
		this.color = color;
		this.radius = getDefaultRadius();
		particles = new ArrayList<Particle>();
	}
	
	
	public int getRadius() {
		return radius;
	}


	public void setRadius(int radius) {
		this.radius = radius;
	}


	public static Family getDefaultgp() {
		return dfamily;
	}
	
	public List<Particle> getParticles() {
		return particles;
	}

	public void addParticle(Particle p)
	{
		particles.add(p);
	}
	
	public void removeParticle(Particle p)
	{
		particles.remove(p);
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public Color getColor() {
		return color;
	}

	public void setColor(Color color) {
		this.color = color;
	}
	
	public static Family getDefaultFamily()
	{
		return dfamily;
	}
	
	public String toString()
	{
		return name;
	}
	
	public static Color[] getSampleColors()
	{
		return colors;
	}
	
	public static Color getColor(String name)
	{
		Color color;
		try {
		    Field field = Class.forName("java.awt.Color").getField(name);
		    color = (Color)field.get(null);//static field, null for parameter
		} catch (Exception e) {
		    color = null; // Not defined
		}
		return color;
	}
	
	public static int getDefaultRadius()
	{
		return 100;
	}
}
