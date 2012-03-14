package xmipp.viewer.models;

import xmipp.jni.MDLabel;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.viewer.ImageDimension;

@SuppressWarnings("serial")
public class MicrographsTable extends MetadataTable {

	public MicrographsTable(GalleryData data) throws Exception {
		super(data);
	}
	
	//Setup columns options
	protected void setupColumns(){
		data.globalRender = true;
		ColumnInfo ctfModel = null;
		for (ColumnInfo ci: data.labels){
			if (ci.allowRender)
				ci.render = true;
			if (ci.getLabel() == MDLabel.MDL_CTFMODEL)
				ctfModel = ci;
			//Remove common prefix CTFCrit_ from columns headers
			ci.labelName = ci.labelName.replace("CTFCrit_", "");
		}
		//Move CTF_MODEL column to the end
		if (ctfModel != null){
			data.labels.remove(ctfModel);
			data.labels.add(ctfModel);
		}
	}
	
	@Override
	protected ImageDimension loadDimension() throws Exception {
		setupColumns();
		return super.loadDimension();
		
	}//function loadDimension

    @Override
    public boolean handleRightClick(int row, int col, XmippPopupMenuCreator xpopup) {
        super.handleRightClick(row, col, xpopup);
        xpopup.setItemVisible(XmippPopupMenuCreator.CTF_PROFILE, true);
        xpopup.setItemVisible(XmippPopupMenuCreator.CTF_RECALCULATE, true);
        return true;
    }
}//class MetadataTable
