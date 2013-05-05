package xmipp.viewer.particlepicker;

import ij.ImagePlus;
import java.awt.Color;
import java.lang.reflect.Field;

import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.ImageGeneric;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.viewer.particlepicker.training.model.TrainingPicker;

public class Family {

	private String name;
	private Color color;
	private int size;
	private int templatesNumber;
	private ImageGeneric templates;
	private String templatesfile;
	private int index;

	
	private static Color[] colors = new Color[] { Color.BLUE, Color.CYAN,
			Color.GREEN, Color.MAGENTA, Color.ORANGE, Color.PINK, Color.YELLOW };
	private static int nextcolor;

	public static Color getNextColor() {
		Color next = colors[nextcolor];
		nextcolor++;
		if (nextcolor == colors.length)
			nextcolor = 0;
		return next;
	}

	public Family(String name, Color color, int size, int templatesNumber, String templatesfile) {
		this(name, color, size, null, templatesNumber, templatesfile);
	}

	public Family(String name, Color color, int size, 
			ParticlePicker ppicker, ImageGeneric templates) {
		if (size < 0 || size > ParticlePicker.fsizemax)
			throw new IllegalArgumentException(String.format(
					"Size should be between 0 and %s, %s not allowed", ParticlePicker.fsizemax,
					size));
		if (name == null || name.equals(""))
			throw new IllegalArgumentException(
					XmippMessage.getEmptyFieldMsg("name"));
		this.name = name;
		this.color = color;
		this.size = size;
		if(templates != null)
			this.templatesfile = templates.getFilename();
		else
			setTemplatesNumber(1);
		if(templates != null)
			try
			{
				templatesNumber = ((int)templates.getNDim());
				this.templates = templates;
				for(int i = 0; i < templatesNumber; i ++)//to initialize templates on c part
					XmippImageConverter.readToImagePlus(templates, ImageGeneric.FIRST_IMAGE + i);
			}
			catch (Exception e)
			{
				e.printStackTrace();
				throw new IllegalArgumentException();
			}
		
		
	}
	

	public Family(String name, Color color, int size, ParticlePicker ppicker, int templatesNumber, String templatesfile) {
		if (size < 0 || size >  ParticlePicker.fsizemax)
			throw new IllegalArgumentException(String.format(
					"Size should be between 0 and %s, %s not allowed",  ParticlePicker.fsizemax,

					size));
		if (name == null || name.equals(""))
			throw new IllegalArgumentException(
					XmippMessage.getEmptyFieldMsg("name"));
		this.name = name;
		this.color = color;
		this.size = size;
		this.templatesfile = templatesfile;
		setTemplatesNumber(templatesNumber);
	
	}
	
	public Family(String name, Color color) {
		this(name, color, getDefaultSize(), null, null);
	}



	public void initTemplates() {

		if(templatesNumber == 0 )
			return;
		try {
			
			this.templates = new ImageGeneric(ImageGeneric.Float);
			templates.resize(size, size, 1, templatesNumber);
			
			templates.write(templatesfile);
			templates.setFilename(templatesfile);
			
			
			
		} catch (Exception e) {
			throw new IllegalArgumentException(e.getMessage());
		}
	}


	public ImageGeneric getTemplates() {
		return templates;
	}
	
	
	
	public ImagePlus getTemplatesImage(long i) {
		try
		{
			ImagePlus imp = XmippImageConverter.convertToImagePlus(templates, i);
			return imp;
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}
	}



	

	public int getSize() {
		return size;
	}

	public void setSize(int size) {
		if (size >  ParticlePicker.fsizemax)
			throw new IllegalArgumentException(String.format(
					"Max size is %s, %s not allowed",  ParticlePicker.fsizemax, size));

		this.size = size;
		initTemplates();
	}


	public String getName() {
		return name;
	}

	public void setName(String name) {
		if (name == null || name.equals(""))
			throw new IllegalArgumentException(
					XmippMessage.getEmptyFieldMsg("name"));
		this.name = name;
	}

	public int getTemplatesNumber() {
		return templatesNumber;
	}

	public void setTemplatesNumber(int num) {
		if(num <= 0)
			throw new IllegalArgumentException(XmippMessage.getIllegalValueMsgWithInfo("Templates Number", Integer.valueOf(num), "Family must have at least one template"));

		this.templatesNumber = num;
		initTemplates();
	}

	public Color getColor() {
		return color;
	}

	public void setColor(Color color) {
		this.color = color;
	}

	public String toString() {
		return name;
	}

	

	public static Color getColor(String name) {
		Color color;
		try {
			Field field = Class.forName("java.awt.Color").getField(name);
			color = (Color) field.get(null);// static field, null for parameter
		} catch (Exception e) {
			color = null; // Not defined
		}
		return color;
	}

	public static int getDefaultSize() {
		return 100;
	}



	public int getRadius() {
		return size / 2;
	}

	public void setTemplate(int index, ImageGeneric ig) {
		this.index = index;
		float[] matrix;
		try {
			//TODO getArrayFloat and setArrayFloat must be call from C both in one function
			matrix = ig.getArrayFloat(ImageGeneric.FIRST_IMAGE,	ImageGeneric.FIRST_SLICE);
			templates.setArrayFloat(matrix, index, ImageGeneric.FIRST_SLICE);
		} catch (Exception e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public String getTemplatesFile()
	{
		return templatesfile;
	}

	public void saveTemplates()
	{
		try
		{
			if(index == templatesNumber)//already filled all initial templates 
				templates.write(getTemplatesFile());
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}
		
	}
	

	

}
