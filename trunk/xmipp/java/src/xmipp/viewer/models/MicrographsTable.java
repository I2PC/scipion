package xmipp.viewer.models;

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
		for (ColumnInfo ci: data.labels)
			if (ci.allowRender)
				ci.render = true;
	}
	
	@Override
	protected ImageDimension loadDimension() throws Exception {
		setupColumns();
		return super.loadDimension();
		
	}//function loadDimension

    @Override
    public boolean handleRightClick(int row, int col, XmippPopupMenuCreator xpopup) {
        super.handleRightClick(row, col, xpopup);
        xpopup.setItemVisible("ShowCTF_mi", true);
        return true;
    }
}//class MetadataTable
