package xmipp.jni;

public class ImageWriteMode {
   public static final int WRITE_READONLY = 0;   //only can read the file /// 
   public static final int WRITE_OVERWRITE = 1; //forget about the old file and overwrite it /// 
   public static final int WRITE_REPLACE = 2;   //replace a particular object by another /// 
   public static final int WRITE_APPEND = 3;    //append and object at the end of a stack = 3; so far can not append stacks /// 
   public static final int WRITE_LAST_LABEL = 4;                       // **** NOTE ****: Do keep this label always at the end /// 
}
