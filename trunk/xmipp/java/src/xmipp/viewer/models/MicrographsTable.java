package xmipp.viewer.models;

import java.util.ArrayList;

import xmipp.jni.MDLabel;
import xmipp.utils.DEBUG;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.viewer.ImageDimension;
import xmipp.viewer.ImageItem;

@SuppressWarnings("serial")
public class MicrographsTable extends MetadataTable {

	protected ArrayList<Integer> busyRows = new ArrayList<Integer>();
	
	public MicrographsTable(GalleryData data) throws Exception {
		super(data);
	}
	
	//Setup columns options
	protected void setupColumns(){
		data.globalRender = true;
		ColumnInfo ctfModel = null;
		ColumnInfo ctfModel2 = null;
		for (ColumnInfo ci: data.labels){
			if (ci.allowRender)
				ci.render = true;
			if (ci.getLabel() == MDLabel.MDL_CTFMODEL)
				ctfModel = ci;
			if (ci.getLabel() == MDLabel.MDL_CTFMODEL2)
				ctfModel2 = ci;
			//Remove common prefix CTFCrit_ from columns headers
			ci.labelName = ci.labelName.replace("CTFCrit_", "");
		}
		//Move CTF_MODEL column to the end
		if (ctfModel != null){
			data.labels.remove(ctfModel);
			data.labels.add(ctfModel);
		}
		if (ctfModel2 != null){
			data.labels.remove(ctfModel2);
			data.labels.add(ctfModel2);
		}
	}
	
	@Override
	protected ImageDimension loadDimension() throws Exception {
		setupColumns();
		return super.loadDimension();
		
	}//function loadDimension
	
	
	@Override
	public Object getValueAt(int row, int col) {
		Object obj = super.getValueAt(row, col);
		
		if (visibleLabels.get(col).render && busyRows.contains(row)) {
			ImageItem ii  = (ImageItem) obj;
			ii.isBusy = true;
		}
		return obj;
	}

    @Override
    public boolean handleRightClick(int row, int col, XmippPopupMenuCreator xpopup) {
        super.handleRightClick(row, col, xpopup);
        xpopup.setItemVisible(XmippPopupMenuCreator.CTF_PROFILE, true);
        xpopup.setItemVisible(XmippPopupMenuCreator.CTF_RECALCULATE, true);
        return true;
    }
    
	public void setRowBusy(int row) {
		busyRows.add(row);
		fireTableRowsUpdated(row, row);
	}

	public void setRowIdle(int row) {
		busyRows.remove(new Integer(row));		
		fireTableRowsUpdated(row, row);
	}
    
    
}//class MetadataTable
