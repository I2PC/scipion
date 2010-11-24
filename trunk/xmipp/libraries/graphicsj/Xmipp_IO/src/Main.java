
import xmipp.io.OpenSPE_;

/*
 * To change this template, choose Tools | Templates
 * and load the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Main {

    public static void main(String args[]) {
        String path = "/home/juanjo/Desktop/imgs_Roberto";
        String name = "img60973.spe";

        OpenSPE_ reader = new OpenSPE_();
        OpenSPE_.load(path, name).show();
        reader.loadThumbnail(path, name, 64, 64).show();
//        reader.loadThumbnail(path, name, 32, 32).show();

        /*        String names[] = Sel_Reader.loadFileNames(path, name);
        String thumbnailFile = names[names.length / 2];
        path = new File(thumbnailFile).getParent();
        name = new File(thumbnailFile).getName();

        Spider_Reader reader = new Spider_Reader();
        reader.load(path, name).show();
        reader.loadThumbnail(path, name, 102, 102).show();
        reader.loadThumbnail(path, name, 204, 204).show();
        reader.loadThumbnail(path, name, 408, 408).show();
        //        reader.loadThumbnail(path, name, 64, 64).show();
        //        reader.loadThumbnail(path, name, 32, 32).show();*/
    }
}
