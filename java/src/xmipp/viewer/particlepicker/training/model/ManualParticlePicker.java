package xmipp.viewer.particlepicker.training.model;

import xmipp.jni.Filename;
import xmipp.jni.MetaData;
import xmipp.viewer.particlepicker.Format;
import xmipp.viewer.particlepicker.IJCommand;
import xmipp.viewer.particlepicker.Micrograph;

	
	public class ManualParticlePicker extends TrainingPicker {
		
		
		public ManualParticlePicker(String selfile, String outputdir, FamilyState mode) {

			this(selfile, outputdir, null, mode);
		}
		
		public ManualParticlePicker(String selfile, String outputdir, String fname, FamilyState mode) {

			super(selfile, outputdir, fname, mode);
			for (TrainingMicrograph m : micrographs)
				loadMicrographData(m);
			setUpdateTemplatesPending(true);
			updateTemplates();
			if(filters.isEmpty())//user just started manual mode and has no filter, I select gaussian blur by default, will be applied when window opens
			{
				filters.add(new IJCommand("Gaussian Blur...", "sigma=2"));
				persistFilters();
			}

		}

		

		public String getImportMicrographName(String path, String filename, Format f)
		{

			String base = Filename.removeExtension(Filename.getBaseName(filename));
			switch (f)
			{
			case Xmipp24:
				return Filename.join(path, base, base + ".raw.Common.pos");
			case Xmipp30:
				return Filename.join(path, base + ".pos");
			case Eman:
				return Filename.join(path, base + ".box");

			default:
				return null;
			}

		}

		
		public Format detectFormat(String path)
		{
			Format[] formats = { Format.Xmipp24, Format.Xmipp30, Format.Eman };
			for (TrainingMicrograph m : micrographs)
			{
				for (Format f : formats)
				{
					if (Filename.exists(getImportMicrographName(path, m.getFile(), f)))
						return f;
				}
			}
			return Format.Unknown;
		}



		/** Return the number of particles imported */
		public int importParticlesFromFolder(String path, Format f, float scale, boolean invertx, boolean inverty)
		{
			if (f == Format.Auto)
				f = detectFormat(path);
			if (f == Format.Unknown)
				return 0;

			String filename;
			int particles = 0;

			for (TrainingMicrograph m : micrographs)
			{
				filename = getImportMicrographName(path, m.getFile(), f);
				System.out.println("  filename: " + filename);
				if (Filename.exists(filename))
					particles += importParticlesFromFile(filename, f, m, scale, invertx, inverty);
			}
			
			return particles;
		}// function importParticlesFromFolder

		

		/** Return the number of particles imported from a file */
		public int importParticlesFromFile(String path, Format f, Micrograph m, float scale, boolean invertx, boolean inverty)
		{
			MetaData md = new MetaData();
			fillParticlesMdFromFile(path, f, m, md, scale, invertx, inverty);
			int particles = (md != null) ? importParticlesFromMd(m, md) : 0;
			saveData(getMicrograph());
			md.destroy();
			return particles;
		}// function importParticlesFromFile
	}