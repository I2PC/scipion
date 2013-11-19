package xmipp.viewer.models;

import java.util.ArrayList;

import xmipp.jni.MDLabel;
import xmipp.jni.MDRow;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.XmippPopupMenuCreator;
import xmipp.viewer.ImageDimension;

@SuppressWarnings("serial")
public class MicrographsTableModel extends MetadataTableModel {

	protected ArrayList<Integer> busyRows = new ArrayList<Integer>();
	
	public MicrographsTableModel(GalleryData data) throws Exception {
		super(data);
	}
	
	//Setup columns options
	protected void setupColumns(){
		data.globalRender = true;
		ColumnInfo ctfModel = null;
		ColumnInfo ctfModel2 = null;
		ColumnInfo firstRender = null;
		for (ColumnInfo ci: data.labels){
			if (ci.allowRender &&
				ci.getLabel() != MDLabel.MDL_MICROGRAPH &&
				ci.getLabel() != MDLabel.MDL_MICROGRAPH_ORIGINAL &&
				ci.getLabel() != MDLabel.MDL_MICROGRAPH_TILTED &&
				ci.getLabel() != MDLabel.MDL_MICROGRAPH_TILTED_ORIGINAL)
			{
				ci.render = true;
				if (firstRender == null)
					firstRender = ci;
			}
			if (ci.getLabel() == MDLabel.MDL_CTF_MODEL)
				ctfModel = ci;
			if (ci.getLabel() == MDLabel.MDL_CTF_MODEL2)
				ctfModel2 = ci;
			//Remove common prefix CTFCrit_ from columns headers
			ci.labelName = ci.labelName.replace("ctfCrit", "");
			
			//Set invisible some colums by default
			if (ci.labelName.startsWith("Psd") || ci.labelName.startsWith("Normality"))
				ci.visible = false;
				
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
		
		//Replace MDL_MICROGRAPH as first render
		if (firstRender != null){
			data.ciFirstRender = firstRender;
			data.zoom = 50;
		}
	}
	
	@Override
	protected ImageDimension loadDimension() throws Exception {
		setupColumns();
		return super.loadDimension();
		
	}//function loadDimension

    @Override
    public boolean handleRightClick(int row, int col, XmippPopupMenuCreator xpopup) {
    	if (isBusy(row, col))
    		return false;
        super.handleRightClick(row, col, xpopup);
        xpopup.setItemVisible(XmippPopupMenuCreator.CTF_PROFILE, true);
        xpopup.setItemVisible(XmippPopupMenuCreator.CTF_RECALCULATE, true);
        return true;
    }
    
	/** Check if the item is busy */
	public boolean isBusy(int row, int col){
		return busyRows.contains(row);
	}
    
	public void setRowBusy(int row) {
		DEBUG.printFormat("setting busy row: %d", row);
		busyRows.add(new Integer(row));
		fireTableRowsUpdated(row, row);
	}

	public void setRowIdle(int row) {
		DEBUG.printFormat("setting idle row: %d", row);
		busyRows.remove(new Integer(row));	
		long objId = data.ids[row];
		String psdFile = data.md.getValueString(MDLabel.MDL_PSD, objId);
		String sortFn = psdFile.replace(".psd", ".tmpSort.xmd");
		MetaData mdRow = new MetaData(sortFn);
		MDRow row2 = new MDRow();
		mdRow.getRow(row2, mdRow.firstObject());
		data.setRow(row2, objId);
		mdRow.destroy();
		refreshRow(row);
	}
	
	/** Function to force the refresh of some row */
	public void refreshRow(int row) {
		try {
			int n = getColumnCount();
			int index;
			String key;
			ColumnInfo ci;
			for (int col = 0; col < n; ++col) {
				ci = data.getColumnInfo(col);
				if (ci.allowRender){
					key = getItemKey(row, ci.getLabel());
					DEBUG.printFormat("Removing item: %s from cache", key);
					cache.remove(key);
				}
			}
			fireTableRowsUpdated(row, row);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
   
    
}//class MetadataTable
