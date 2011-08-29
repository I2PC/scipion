package model;
import java.awt.Color;
import java.io.File;
import java.lang.reflect.Field;



public class Family {
	
	private String name;
	private Color color;
	private int size;
	private FamilyState step;
	int particles = 0;
	
	public int getParticlesNumber()
	{
		return particles;
	}


	private static int sizemax = 1000;
	private static Family dfamily = new Family("Default", Color.green);
	private static Color[] colors = new Color[]{Color.BLUE, Color.CYAN, 
										Color.GREEN,
										Color.MAGENTA, Color.ORANGE, 
										Color.PINK, Color.RED, Color.YELLOW};
	private static int nextcolor;
	
	
	
	
	
	public static Color getNextColor()
	{
		Color next = colors[nextcolor];
		nextcolor ++;
		if(nextcolor == colors.length)
			nextcolor = 0;
		return next;
	}
	
	
	public Family(String name, Color color, int size)
	{
		this(name, color, size, FamilyState.Manual);
	}
	
	public Family(String name, Color color, int size, FamilyState state)
	{
		if(size < 0 || size > sizemax)
			throw new IllegalArgumentException(String.format("Size should be between 0 and %s, %s not allowed", sizemax, size));
		if (name == null || name.equals(""))
			throw new IllegalArgumentException(Constants.getEmptyFieldMsg("name"));
		this.name = name;
		this.color = color;
		this.size = size;
		this.step = state;
	}
	
	public String getOTrainingFilename()
	{
		return ParticlePicker.getOutputPath(String.format("%s_%s", name, ParticlePicker.getTrainingFilenameGeneric()));
	}
	
	public String getOTrainingMaskFilename()
	{
		return ParticlePicker.getOutputPath(String.format("%s_%s", name, ParticlePicker.getTrainingMaskFilenameGeneric()));
	}
	
	
	public Family(String name, Color color)
	{
		this(name, color, getDefaultSize(), FamilyState.Manual);
	}
	
	
	public FamilyState getStep()
	{
		return step;
	}
	
	public String getOutputRoot()
	{
		return ParticlePicker.getOutputPath(name);
	}
	
	public void goToNextStep()
	{
		validateNextStep();
		this.step = ParticlePicker.nextStep(step);
		ParticlePicker.getInstance().resetModel(this);
		
	}
	
	
	
	public void validateNextStep()
	{
		int min = ParticlePicker.getMinForTraining();
		FamilyState next = ParticlePicker.nextStep(step);
		if(next == FamilyState.Supervised && particles < min)
			throw new IllegalArgumentException(String.format("You should have at least %s particles to go to %s mode", min, FamilyState.Supervised));
		if(!ParticlePicker.getInstance().hasEmptyMicrographs(this))
			throw new IllegalArgumentException(String.format("There are no available micrographs for %s step", FamilyState.Supervised));
		
	}
	
	public static String getOFilename()
	{
		return ParticlePicker.getOutputPath("families.xmd");
	}
	
	public int getSize() {
		return size;
	}


	public void setSize(int size) {
		if(size > sizemax)
			throw new IllegalArgumentException(String.format("Max size is %s, %s not allowed", sizemax, size));
		this.size = size;
	}


	public static Family getDefaultgp() {
		return dfamily;
	}
	

	public String getName() {
		return name;
	}

	public void setName(String name) {
		if (name == null || name.equals(""))
			throw new IllegalArgumentException(Constants.getEmptyFieldMsg("name"));
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
	
	public static int getDefaultSize()
	{
		return 100;
	}
}
