/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.models;

import javax.swing.table.DefaultTableModel;

import ij.IJ;
import xmipp.utils.Cache;
import xmipp.utils.DEBUG;
// import xmipp.viewer.imageitems.tableitems.GalleryImageItem;
import ij.ImagePlus;
// import xmipp.viewer.metadata.images.TableFileItem;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class SimpleMetadataTableModel extends DefaultTableModel {

    private int ENABLED_COLUMN_INDEX;
    private int MD_LABELS[];
    private String TEXT_LABELS[];
    private String filename;
    private String blocks[];
    private int selectedBlock = 0;
    private long ids[];
    protected static Cache<String, ImagePlus> cache = new Cache<String, ImagePlus>();
    //protected boolean containsImageLabel = false;
    protected boolean containsImages = false;

    public SimpleMetadataTableModel(String filename) {
        super();
        if(filename != null){
	        this.filename = filename;
	
	        loadBlocks(filename);
	        selectBlock(Filename.getPrefix(filename));
        }
    }

    private void selectBlock(String selectedBlock) {
        int n = 0;

        if (selectedBlock != null) {
            n = findBlock(selectedBlock, blocks);
        }

        selectBlock(n > 0 ? n : 0);
    }

    public void selectBlock(int selectedBlock) {
        this.selectedBlock = selectedBlock;
    }

    public int getSelectedBlockIndex() {
        return selectedBlock;
    }

    public String getSelectedBlock() {
        return blocks[getSelectedBlockIndex()];
    }

//    public boolean containsImageLabel() {
//        return containsImageLabel;
//    }
    public boolean containsImages() {
        return containsImages;
    }

    public String getFilename() {
        return Filename.getFilename(filename);
    }

    public String getPath() {
        String block = getSelectedBlock();

        return (!block.isEmpty() ? block + Filename.SEPARATOR : "") + getFilename();
    }

    private void clear() {
        cache.clear();

        getDataVector().removeAllElements();

        setColumnCount(0);
        ENABLED_COLUMN_INDEX = -1;
    }

    public String[] getBlocks() {
        return blocks;
    }

    private void loadBlocks(String filename) {
        try {
            String name = Filename.getFilename(filename);
            String block = Filename.getPrefix(filename);

            blocks = MetaData.getBlocksInMetaDataFile(name);

            // No blocks, so set at least one item for the dropdown list.
            if (blocks.length < 1) {
                blocks = new String[]{""};
            }

            if (block != null) {
                // Check if selected block is in the list, otherwise adds it (like in 'class_[regexp]')
                int n = findBlock(block, blocks);
                if (n < 0) {
                    String blocks_[] = new String[blocks.length + 1];
                    blocks_[0] = block;
                    System.arraycopy(blocks, 0, blocks_, 1, blocks.length);
                    blocks = blocks_;
                }
            }
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }
    }

    public static int findBlock(String block, String blocks[]) {
        for (int n = 0; n < blocks.length; n++) {
            if (blocks[n].equals(block)) {
                return n;
            }
        }
        return -1;
    }

    public void reload() {
        load(getPath());
    }

    public void initFromMetadata(MetaData md) throws Exception{
        // Contains field enabled ?
        boolean hasEnabledField = true;
        if (!md.containsLabel(MDLabel.MDL_ENABLED)) {
            if (md.containsLabel(MDLabel.MDL_IMAGE)) {
                md.addLabel(MDLabel.MDL_ENABLED);
            }
            hasEnabledField = false;
        }

        // Contains image (for gallery)?
        if (md.containsLabel(MDLabel.MDL_IMAGE)) {
            containsImages = true;

            setCacheSize(md);
        }

        // Store ids.
        ids = md.findObjects();

        // Reads all rows.
        MD_LABELS = md.getActiveLabels();
        // Builds columns structure.
        TEXT_LABELS = new String[MD_LABELS.length];
        for (int i = 0; i < MD_LABELS.length; i++) {
            TEXT_LABELS[i] = MetaData.label2Str(MD_LABELS[i]);
            addColumn(TEXT_LABELS[i]);
        }
        Object row[] = new Object[MD_LABELS.length];

        for (long id : ids) {
            // Field enabled does not exist, so adds it (true by default).
            if (!hasEnabledField) {
                md.setValueInt(MDLabel.MDL_ENABLED, 1, id);
            }

            for (int i = 0; i < MD_LABELS.length; i++) {
                int label = MD_LABELS[i];
                Class class_ = MetaData.getLabelClass(label);

                if (class_ == String.class) {
                    String value = md.getValueString(label, id);
                    row[i] = value;
                    if (MetaData.isImage(label)) {
                        String path = md.getValueString(label, id, true);
                        //row[i] = new GalleryImageItem(path, value, cache);
                        //row[i] = new GalleryImageItem(id, md, label, cache);
                        containsImages = true;
                    } else if (MetaData.isMetadata(label)) {
                        String path_ = md.getValueString(label, id, true);
                        //row[i] = new TableFileItem(path_, value);
                        row[i] = path_;
                    } else if (MetaData.isTextFile(label)) {
                        String path_ = md.getValueString(label, id, true);
                        //row[i] = new TableFileItem(path_, value);
                        row[i] = path_;
                    }
                } else if (class_ == Double.class) {
                    row[i] = md.getValueDouble(label, id);
                } else if (class_ == Long.class) {
                    row[i] = md.getValueLong(label, id);
                } else if (class_ == Integer.class) {
                    int value = md.getValueInt(label, id);

                    // Special case.
                    if (label == MDLabel.MDL_ENABLED) {
                        row[i] = value == 1;
                        ENABLED_COLUMN_INDEX = i;   // Store it for future changes.
                    } else {
                        row[i] = value;
                    }
                } else if (class_ == Boolean.class) {
                    row[i] = md.getValueBoolean(label, id);
                } else {
                    row[i] = "Unknown type";
                }
            }

            // Adds a new row data to table.
            addRow(row);
        }
    }
    
    
    public void load(String filename) {
        try {
            clear();    // Clear the whole data.

            MetaData md = new MetaData(filename);
            initFromMetadata(md);
            md.destroy();

        } catch (Exception ex) {
            DEBUG.printException(ex);
            IJ.error(ex.getMessage());
        }
    }

    void setCacheSize(MetaData md) throws Exception {
        // Calculates cache elements size.
        String firstImage = md.getValueString(MDLabel.MDL_IMAGE, md.firstObject(), true);
        ImageGeneric image = new ImageGeneric(firstImage);

        int imageSize = image.getXDim() * image.getYDim() * Cache.MAXPXSIZE;
        int elements = Cache.MEMORY_SIZE / imageSize;

        cache.resize(elements > 0 ? elements : 1);
    }

    public String[] getLabels() {
        return TEXT_LABELS;
    }

    public int getEnabledColumnIndex() {
        return ENABLED_COLUMN_INDEX;
    }

    @Override
    public Class getColumnClass(int column) {
        Object item = null;
        try {
            item = getValueAt(0, column);
        } catch (Exception e) {
        }

        return item != null ? item.getClass() : Object.class;
    }

    @Override
    public int getColumnCount() {
        return MD_LABELS != null ? MD_LABELS.length : 0;
    }

    public boolean isRowEnabled(int row) {
        Object value = ENABLED_COLUMN_INDEX >= 0 ? getValueAt(row, ENABLED_COLUMN_INDEX) : null;

        return value == null || ((Boolean) value).booleanValue();
    }

    public void enableAllRows(boolean enabled) {
        // Updates table.
        for (int i = 0; i < getRowCount(); i++) {
            //setRowEnabled(i, enabled);
            setValueAt(enabled, i, ENABLED_COLUMN_INDEX);
        }
    }

    public boolean save(String fileName) {
    	/*
        try {
            MetaData md = new MetaData();

            for (int i = 0; i < MD_LABELS.length; i++) {
                md.addLabel(MD_LABELS[i]);
            }

            for (int row = 0; row < getRowCount(); row++) {
                long id = md.addObject();

                for (int column = 0; column < getColumnCount(); column++) {
                    int label = MD_LABELS[column];
                    Class class_ = getColumnClass(column);
                    Object item = getValueAt(row, column);

                    if (class_ == TableFileItem.class) {
                        md.setValueString(label, ((TableFileItem) item).getOriginalValue(), id);
                        //} else if (class_ == GalleryImageItem.class) {
                    } else if (class_ == GalleryImageItem.class) {
                        md.setValueString(label, ((GalleryImageItem) item).getOriginalValue(), id);
                    } else if (class_ == String.class) {
                        md.setValueString(label, ((String) item), id);
                    } else if (class_ == Double.class) {
                        md.setValueDouble(label, ((Double) item).doubleValue(), id);
                    } else if (class_ == Boolean.class) {
                        if (column == ENABLED_COLUMN_INDEX) {
                            md.setValueInt(label, ((Boolean) item) == Boolean.TRUE ? 1 : 0, id);
                        } else {
                            md.setValueBoolean(label, ((Boolean) item).booleanValue(), id);
                        }
                    } else if (class_ == Integer.class) {
                        md.setValueInt(label, ((Integer) item).intValue(), id);
                    }
                }
            }

            md.write(fileName);

            return true;
        } catch (Exception ex) {
            DEBUG.printException(ex);
            IJ.error(ex.getMessage());
        }

        return true;*/
    	return false;
    }
    	
    
}
