package trainingpicker.model;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import xmipp.Particle;
import xmipp.MDLabel;
import xmipp.MetaData;


public class ReviewParticlePicker extends TrainingPicker {
	
	private String reviewfile;
	private Family reviewfamily;
	
	public String getReviewFile() {
		return reviewfile;
	}
	
	public Family getReviewFamily() {
		return reviewfamily;
	}
	

	public ReviewParticlePicker(String selfile, String outputdir, String reviewfile) 
	{
		super(selfile, outputdir, FamilyState.Review);
		if (!new File(reviewfile).exists())
			throw new IllegalArgumentException(
					Constants.getNoSuchFieldValueMsg("review file", reviewfile));
		this.reviewfile = reviewfile;
		String[] parts = reviewfile.split(File.separator);
		String familyname = parts[parts.length - 1].split("_")[0];
		this.reviewfamily = getFamily(familyname);
		if (reviewfamily == null)
			throw new IllegalArgumentException(
					Constants.getNoSuchFieldValueMsg("family", familyname));
		
		loadMicrographs();
		families.clear();
		families.add(reviewfamily);
	}
	
	
	
	
	@Override
	public void persistMicrographs() {
		exportData(reviewfamily);
	}
	
	

	@Override
	public void loadMicrographs() {
		micrographs.clear();
		TrainingMicrograph micrograph;
		String ctf = null, filename;
		try {
			MetaData md = new MetaData(getMicrographsSelFile());
			boolean existsctf = md.containsLabel(MDLabel.MDL_PSD_ENHANCED);
			long[] ids = md.findObjects();
			for (long id : ids) {

				filename = md.getValueString(MDLabel.MDL_IMAGE, id);
				if (existsctf)
					ctf = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				micrograph = new TrainingMicrograph(filename, ctf, families, getMode());
				micrographs.add(micrograph);
			}
			if (micrographs.size() == 0)
				throw new IllegalArgumentException(String.format(
						"No micrographs specified on %s", getMicrographsSelFile()));
			int x, y;
			double cost;
			TrainingParticle particle;
			List<String> blocks = Arrays.asList(MetaData
					.getBlocksInMetaDataFile(reviewfile));
			String block;
			for (TrainingMicrograph m : micrographs) {
				block = "mic_" + m.getName();
				if (blocks.contains(block)) {
					md = new MetaData(block + "@" + reviewfile);

					ids = md.findObjects();
					for (long id : ids) {

						x = md.getValueInt(MDLabel.MDL_XINT, id);
						y = md.getValueInt(MDLabel.MDL_YINT, id);
						cost = md.getValueDouble(MDLabel.MDL_COST, id);
						if(cost > 1)
							m.addManualParticle(new TrainingParticle(x, y, reviewfamily, m, cost));
						else
							m.addAutomaticParticle(new AutomaticParticle(x, y, reviewfamily, m, cost, false));
					}
				}
			}

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}
	



}
