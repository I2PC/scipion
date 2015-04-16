package xmipp.jni;

public class CastWriteMode {
   public static final int CW_CAST = 0;       //Only cast the data type /// 
   public static final int CW_CONVERT = 1;    //Convert the data from one type to another /// 
   public static final int CW_ADJUST = 2;     //Adjust the histogram to fill the gray level range /// 
   public static final int CW_LAST_LABEL = 3;                       // **** NOTE ****: Do keep this label always at the end /// 
}
