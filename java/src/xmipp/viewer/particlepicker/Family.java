package xmipp.viewer.particlepicker;

import ij.ImagePlus;
import ij.ImageStack;

import java.awt.Color;
import java.lang.reflect.Field;

import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.ImageGeneric;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.training.model.FamilyState;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.viewer.particlepicker.training.model.TrainingPicker;

public class Family {

	private String name;
	private Color color;
	private int size;
	private FamilyState state;
	private int templatesNumber;
	private ImageGeneric templates;
	private String templatesfile;
	protected boolean updateTemplatesPending;

	
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
		this(name, color, size, FamilyState.Manual, null, templatesNumber, templatesfile);
	}

	public Family(String name, Color color, int size, FamilyState state,
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
		this.state = state;
		this.templatesfile = templates.getFilename();
		if(templates == null)
			setTemplatesNumber(1);
		else
			try
			{
				templatesNumber = ((int)templates.getNDim());
				this.templates = templates;
				for(int i = 0; i < templatesNumber; i ++)//to initialize templates on c part
					getTemplatesImage(ImageGeneric.FIRST_IMAGE + i);
			}
			catch (Exception e)
			{
				e.printStackTrace();
				throw new IllegalArgumentException();
			}
		
		
	}
	

	public Family(String name, Color color, int size, FamilyState state,
			ParticlePicker ppicker, int templatesNumber, String templatesfile) {
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
		this.state = state;
		this.templatesfile = templatesfile;
		setTemplatesNumber(templatesNumber);
	
	}
	
	public String getTemplatesFile()
	{
		return templatesfile;
	}
	
	public void initTemplates()
	{
		if(templatesNumber == 0)
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
			ImagePlus imp = XmippImageConverter.readToImagePlus(templates, i);
			return imp;
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}
	}


	public FamilyState getStep() {
		return state;
	}

	public void goToNextStep(TrainingPicker ppicker) {
		validateNextStep(ppicker);
		this.state = TrainingPicker.nextStep(state);
	}

	public void goToPreviousStep() {
		this.state = TrainingPicker.previousStep(state);
	}

	public void validateNextStep(TrainingPicker ppicker) {
		int min = SupervisedParticlePicker.getMinForTraining();
		FamilyState next = TrainingPicker.nextStep(state);
		if (next == FamilyState.Supervised
				&& ppicker.getManualParticlesNumber(this) < min)
			throw new IllegalArgumentException(String.format(
					"You should have at least %s particles to go to %s mode",
					min, FamilyState.Supervised));
		if (!ppicker.hasEmptyMicrographs(this) && next != FamilyState.Review)
			throw new IllegalArgumentException(String.format(
					"There are no available micrographs for %s step",
					FamilyState.Supervised));

	}

	// public static String getOFilename()
	// {
	// return ParticlePicker.getInstance().getOutputPath("families.xmd");
	// }

	public int getSize() {
		return size;
	}

	public void setSize(int size) {
		if (size > ParticlePicker.fsizemax)
			throw new IllegalArgumentException(String.format(
					"Max size is %s, %s not allowed", ParticlePicker.fsizemax, size));
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
		setUpdateTemplatesPending(true);
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

	public static Color[] getSampleColors() {
		return colors;
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

	public void setState(FamilyState state) {
		this.state = state;

	}

	public int getRadius() {
		return size / 2;
	}

	public void setTemplate(int index, ImageGeneric ig) {
		
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
	
	public void setUpdateTemplatesPending(boolean b)
	{
		updateTemplatesPending = b;

	}

	public boolean getUpdateTemplatesPending()
	{
		return updateTemplatesPending;
	}

}
