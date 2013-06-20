package xmipp.viewer.particlepicker;

import ij.ImagePlus;
import java.awt.Color;
import java.io.File;
import java.lang.reflect.Field;
import java.util.List;

import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.ImageGeneric;
import xmipp.jni.Particle;
import xmipp.utils.TasksManager;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.training.model.FamilyState;
import xmipp.viewer.particlepicker.training.model.MicrographFamilyData;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.viewer.particlepicker.training.model.TrainingMicrograph;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;
import xmipp.viewer.particlepicker.training.model.TrainingPicker;

public class Family
{

	private String name;
	private Color color;
	private int size;
	private FamilyState state;
	private int templatesNumber;
	private ImageGeneric templates;
	private String templatesfile;
	private int templateindex;
	private ParticlePicker picker;

	private static Color[] colors = new Color[] { Color.BLUE, Color.CYAN, Color.GREEN, Color.MAGENTA, Color.ORANGE, Color.PINK, Color.YELLOW };
	private static int nextcolor;

	public static Color getNextColor()
	{
		Color next = colors[nextcolor];
		nextcolor++;
		if (nextcolor == colors.length)
			nextcolor = 0;
		return next;
	}

	public Family(String name, Color color, int size, ParticlePicker ppicker, int templatesNumber)
	{
		this(name, color, size, FamilyState.Manual, ppicker, templatesNumber);
	}
	

	public Family(String name, Color color, int size, FamilyState state, ParticlePicker ppicker, int templatesNumber)
	{
		if (size < 0 || size > ParticlePicker.fsizemax)
			throw new IllegalArgumentException(String.format("Size should be between 0 and %s, %s not allowed", ParticlePicker.fsizemax, size));
		if (name == null || name.equals(""))
			throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("name"));
		try
		{
			this.picker = ppicker;
			this.name = name;
			this.color = color;
			this.size = size;
			this.state = state;
			this.templatesNumber = templatesNumber;
			this.templatesfile = picker.getTemplatesFile(name);
			if (new File(templatesfile).exists())
			{

				this.templates = new ImageGeneric(templatesfile);
				templates.read(templatesfile, false);
				templateindex = templatesNumber;//all images read

			}
			else
				initTemplates();
			


			
		}
		catch (Exception e)
		{
			e.printStackTrace();
			throw new IllegalArgumentException();
		}
	}


	public String getTemplatesFile()
	{
		return templatesfile;
	}

	public synchronized void initTemplates()
	{
		if (templatesNumber == 0)
			return;
		try
		{
			this.templates = new ImageGeneric(ImageGeneric.Float);
			templates.resize(size, size, 1, templatesNumber);
			templates.write(templatesfile);
			templateindex = 0;
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}
	}
	
	public void updateTemplates(TrainingPicker picker)
	{
		if (getStep() != FamilyState.Manual)
			return;

		initTemplates();
		ImageGeneric igp;
		List<TrainingParticle> particles;
		MicrographFamilyData mfd;
		TrainingParticle particle;
		double[] align;
		try
		{
			for (TrainingMicrograph m : picker.getMicrographs())
			{
				mfd = m.getFamilyData(this);
				for (int i = 0; i < mfd.getManualParticles().size(); i++)
				{
					particles = mfd.getManualParticles();
					particle = particles.get(i);
					igp = particle.getImageGeneric();
					if (getTemplateIndex() < getTemplatesNumber())
						setTemplate(igp);
					else
					{
						align = getTemplates().alignImage(igp);
						particle.setLastalign(align);
					}
				}
			}
			saveTemplates();
			
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}
	}


	public ImageGeneric getTemplates()
	{
		return templates;
	}



	public FamilyState getStep()
	{
		return state;
	}

	public void goToNextStep(TrainingPicker ppicker)
	{
		validateNextStep(ppicker);
		this.state = TrainingPicker.nextStep(state);
	}

	public void goToPreviousStep()
	{
		this.state = TrainingPicker.previousStep(state);
	}

	public void validateNextStep(TrainingPicker ppicker)
	{
		int min = SupervisedParticlePicker.getMinForTraining();
		FamilyState next = TrainingPicker.nextStep(state);
		if (next == FamilyState.Supervised && ppicker.getManualParticlesNumber(this) < min)
			throw new IllegalArgumentException(String.format("You should have at least %s particles to go to %s mode", min, FamilyState.Supervised));
		if (!ppicker.hasEmptyMicrographs(this) && next != FamilyState.Review)
			throw new IllegalArgumentException(String.format("There are no available micrographs for %s step", FamilyState.Supervised));

	}

	public int getSize()
	{
		return size;
	}

	public void setSize(int size)
	{
		if (size > ParticlePicker.fsizemax)
			throw new IllegalArgumentException(String.format("Max size is %s, %s not allowed", ParticlePicker.fsizemax, size));
		this.size = size;

	}


	public String getName()
	{
		return name;
	}

	public void setName(String name)
	{
		if (name == null || name.equals(""))
			throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("name"));
		this.name = name;
	}

	public int getTemplatesNumber()
	{
		return templatesNumber;
	}

	public void setTemplatesNumber(int num) {
		if(num <= 0)
			throw new IllegalArgumentException(XmippMessage.getIllegalValueMsgWithInfo("Templates Number", Integer.valueOf(num), "Family must have at least one template"));

		this.templatesNumber = num;

	}

	public Color getColor()
	{
		return color;
	}

	public void setColor(Color color)
	{
		this.color = color;
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
		try
		{
			Field field = Class.forName("java.awt.Color").getField(name);
			color = (Color) field.get(null);// static field, null for parameter
		}
		catch (Exception e)
		{
			color = null; // Not defined
		}
		return color;
	}

	public static int getDefaultSize()
	{
		return 100;
	}

	public void setState(FamilyState state)
	{
		this.state = state;

	}

	public int getRadius()
	{
		return size / 2;
	}

	public synchronized void setTemplate(ImageGeneric ig)
	{

		float[] matrix;
		try
		{
			// TODO getArrayFloat and setArrayFloat must be call from C both in
			// one function
			matrix = ig.getArrayFloat(ImageGeneric.FIRST_IMAGE, ImageGeneric.ALL_SLICES);
			templates.setArrayFloat(matrix, ImageGeneric.FIRST_IMAGE + templateindex, ImageGeneric.ALL_SLICES);
			templateindex++;

		}
		catch (Exception e)
		{
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
	}


	public synchronized void saveTemplates()
	{
		try
		{
			templates.write(templatesfile);
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}
	}



	public int getTemplateIndex()
	{
		return templateindex;
	}

	public synchronized void centerParticle(TrainingParticle p)
	{

		if (((TrainingPicker)picker).getManualParticlesNumber(this) < templatesNumber)
			return;//missing templates
		Particle shift = null;
		try
		{
			ImageGeneric igp = p.getImageGeneric();
			shift = templates.bestShift(igp);
			double distance = Math.sqrt(Math.pow(shift.getX(),  2) + Math.pow(shift.getY(), 2))/size;
			System.out.printf("normalized distance:%.2f\n", distance);
			if(distance < 0.25)
			{
				p.setX(p.getX() + shift.getX());
				p.setY(p.getY() + shift.getY());
			}
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public synchronized void applyAlignment(TrainingParticle particle, ImageGeneric igp, double[] align)
	{
		try
		{
			particle.setLastalign(align);
			templates.applyAlignment(igp, particle.getTemplateIndex(), particle.getTemplateRotation(), 
					particle.getTemplateTilt(), 
					particle.getTemplatePsi());
			//System.out.printf("adding particle: %d %.2f %.2f %.2f\n", particle.getTemplateIndex(), particle.getTemplateRotation(), particle.getTemplateTilt(), particle.getTemplatePsi());
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	

	

}
