package xmipp.particlepicker;

import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.io.File;
import java.lang.reflect.Field;

import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.ImageGeneric;
import xmipp.particlepicker.training.model.FamilyState;
import xmipp.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.particlepicker.training.model.TrainingPicker;
import xmipp.utils.XmippMessage;

public class Family {

	private String name;
	private Color color;
	private int size;
	private FamilyState state;
	private int templatesNumber;
	private ImageGeneric templates;

	private static int sizemax = 1000;
	private static Family dfamily = new Family("DefaultFamily", Color.green, sizemax/2, 1);
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

	public Family(String name, Color color, int size, int templatesNumber) {
		this(name, color, size, FamilyState.Manual, null, templatesNumber);
	}

	public Family(String name, Color color, int size, FamilyState state,
			ParticlePicker ppicker, ImageGeneric templates) {
		if (size < 0 || size > sizemax)
			throw new IllegalArgumentException(String.format(
					"Size should be between 0 and %s, %s not allowed", sizemax,
					size));
		if (name == null || name.equals(""))
			throw new IllegalArgumentException(
					XmippMessage.getEmptyFieldMsg("name"));
		this.name = name;
		this.color = color;
		this.size = size;
		this.state = state;
		if(templates != null)
			try
			{
				templatesNumber = ((int)templates.getNDim());
				this.templates = templates;
				XmippImageConverter.readToImagePlus(templates);//temporarily to avoid error
			}
			catch (Exception e)
			{
				e.printStackTrace();
				throw new IllegalArgumentException();
			}
		else
			setTemplatesNumber(1);
		
	}
	

	public Family(String name, Color color, int size, FamilyState state,
			ParticlePicker ppicker, int templatesNumber) {
		if (size < 0 || size > sizemax)
			throw new IllegalArgumentException(String.format(
					"Size should be between 0 and %s, %s not allowed", sizemax,
					size));
		if (name == null || name.equals(""))
			throw new IllegalArgumentException(
					XmippMessage.getEmptyFieldMsg("name"));
		this.name = name;
		this.color = color;
		this.size = size;
		this.state = state;
		setTemplatesNumber(templatesNumber);
	
	}
	
	public void initTemplates()
	{
		if(templatesNumber == 0)
			return;
		try {
			this.templates = new ImageGeneric(ImageGeneric.Float);
			templates.resize(size, size, 1, templatesNumber);
			System.out.println("NDim on initTemplates: " + templates.getNDim());
		} catch (Exception e) {
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public ImageGeneric getTemplates() {
		return templates;
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
		if (size > sizemax)
			throw new IllegalArgumentException(String.format(
					"Max size is %s, %s not allowed", sizemax, size));
		this.size = size;
		initTemplates();
	}

	public static Family getDefaultgp() {
		return dfamily;
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

	public static Family getDefaultFamily() {
		return dfamily;
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
			ig.printShape();
			matrix = ig.getArrayFloat(ImageGeneric.FIRST_IMAGE,
					ImageGeneric.FIRST_SLICE);
			templates.setArrayFloat(matrix, index, ImageGeneric.FIRST_SLICE);
//			templates.printShape();
		} catch (Exception e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
	}
}
