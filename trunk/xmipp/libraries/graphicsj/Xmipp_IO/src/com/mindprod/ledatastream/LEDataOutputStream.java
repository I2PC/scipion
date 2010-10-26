/*
 * @(#)LEDataOutputStream.java
 *
 * Summary: Little-endian version of DataOutputStream.
 *
 * Copyright: (c) 1998-2009 Roedy Green, Canadian Mind Products, http://mindprod.com
 *
 * Licence: This software may be copied and used freely for any purpose but military.
 *          see http://mindprod.com/contact/nonmil.html
 *
 * Requires: JDK 1.1+
 *
 * Created with: IntelliJ IDEA IDE.
 *
 * Version History:
 *  1.0 1998-01-06
 *  1.1 1998-01-07 - officially implements DataInput
 *  1.2 1998-01-09 - add LERandomAccessFile
 *  1.3 1998-08-28
 *  1.4 1998-11-10 - add new address and phone.
 *  1.5 1999-10-08 - use com.mindprod.ledatastream package name.
 *  1.6 2005-06-13 - made readLine deprecated
 *  1.7 2007-01-01
 *  1.8 2007-05-24 - add pad, icon, pass Intellij inspector
 */

package com.mindprod.ledatastream;

import java.io.DataOutput;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/**
 * Little-endian version of DataOutputStream.
 *
 * Very similar to
 * DataOutputStream except it writes little-endian instead of
 * big-endian binary data. We can't extend DataOutputStream
 * directly since it has only final methods. This forces us
 * implement LEDataOutputStream with a DataOutputStream object,
 * and use wrapper methods.
 *
 * @author Roedy Green, Canadian Mind Products
 * @since 1998-01-06
 * @version 1.8 2007-05-24
 */

public final class LEDataOutputStream implements DataOutput
    {
    // ------------------------------ CONSTANTS ------------------------------

    /**
     * undisplayed copyright notice.
     *
     * @noinspection UnusedDeclaration
     */
    private static final String EMBEDDED_COPYRIGHT =
            "copyright (c) 1999-2009 Roedy Green, Canadian Mind Products, http://mindprod.com";

    // ------------------------------ FIELDS ------------------------------

    /**
     * to get at big-Endian write methods of DataOutPutStream.
     *
     * @noinspection WeakerAccess
     */
    protected final DataOutputStream dis;

    /**
     * work array for composing output.
     *
     * @noinspection WeakerAccess
     */
    protected final byte[] work;

    // -------------------------- PUBLIC INSTANCE  METHODS --------------------------

    /**
     * constructor.
     *
     * @param out the outputstream we write little endian binary data onto.
     */
    public LEDataOutputStream( OutputStream out )
        {
        this.dis = new DataOutputStream( out );
        work = new byte[8];// work array for composing output
        }

    /**
     * Close stream.
     *
     * @throws IOException if close fails.
     */
    public final void close() throws IOException
        {
        dis.close();
        }

    /**
     * Flush stream without closing.
     *
     * @throws IOException if flush fails.
     */
    public void flush() throws IOException
        {
        dis.flush();
        }

    /**
     * Get size of stream.
     *
     * @return bytes written so far in the stream. Note this is a int, not a long as you would exect. This because the
     *         underlying DataInputStream has a design flaw.
     */
    public final int size()
        {
        return dis.size();
        }

    /**
     * This method writes only one byte, even though it says int (non-Javadoc)
     *
     * @param ib the byte to write.
     *
     * @throws IOException if write fails.
     * @see java.io.DataOutput#write(int)
     */
    public final synchronized void write( int ib ) throws IOException
        {
        dis.write( ib );
        }

    /**
     * Write out an array of bytes.
     *
     * @throws IOException if write fails.
     * @see java.io.DataOutput#write(byte[])
     */
    public final void write( byte ba[] ) throws IOException
        {
        dis.write( ba, 0, ba.length );
        }

    /**
     * Writes out part of an array of bytes.
     *
     * @throws IOException if write fails.
     * @see java.io.DataOutput#write(byte[],int,int)
     */
    public final synchronized void write( byte ba[],
                                          int off,
                                          int len ) throws IOException
        {
        dis.write( ba, off, len );
        }

    /**
     * Write a booleans as one byte.
     *
     * @param v boolean to write.
     *
     * @throws IOException if write fails.
     * @see java.io.DataOutput#writeBoolean(boolean)
     */
    /* Only writes one byte */
    public final void writeBoolean( boolean v ) throws IOException
        {
        dis.writeBoolean( v );
        }

    /**
     * write a byte.
     *
     * @param v the byte to write.
     *
     * @throws IOException if write fails.
     * @see java.io.DataOutput#writeByte(int)
     */
    public final void writeByte( int v ) throws IOException
        {
        dis.writeByte( v );
        }

    /**
     * Write a string.
     *
     * @param s the string to write.
     *
     * @throws IOException if write fails.
     * @see java.io.DataOutput#writeBytes(java.lang.String)
     */
    public final void writeBytes( String s ) throws IOException
        {
        dis.writeBytes( s );
        }

    /**
     * Write a char. Like DataOutputStream.writeChar. Note the parm is an int even though this as a writeChar
     *
     * @param v the char to write
     *
     * @throws IOException if write fails.
     */
    public final void writeChar( int v ) throws IOException
        {
        // same code as writeShort
        work[ 0 ] = ( byte ) v;
        work[ 1 ] = ( byte ) ( v >> 8 );
        dis.write( work, 0, 2 );
        }

    /**
     * Write a string, not a char[]. Like DataOutputStream.writeChars, flip endianness of each char.
     *
     * @throws IOException if write fails.
     */
    public final void writeChars( String s ) throws IOException
        {
        int len = s.length();
        for ( int i = 0; i < len; i++ )
            {
            writeChar( s.charAt( i ) );
            }
        }// end writeChars

    /**
     * Write a double.
     *
     * @param v the double to write. Like DataOutputStream.writeDouble.
     *
     * @throws IOException if write fails.
     */
    public final void writeDouble( double v ) throws IOException
        {
        writeLong( Double.doubleToLongBits( v ) );
        }

    /**
     * Write a float. Like DataOutputStream.writeFloat.
     *
     * @param v the float to write.
     *
     * @throws IOException if write fails.
     */
    public final void writeFloat( float v ) throws IOException
        {
        writeInt( Float.floatToIntBits( v ) );
        }

    /**
     * Write an int, 32-bits.  Like DataOutputStream.writeInt.
     *
     * @param v the int to write
     *
     * @throws IOException if write fails.
     */
    public final void writeInt( int v ) throws IOException
        {
        work[ 0 ] = ( byte ) v;
        work[ 1 ] = ( byte ) ( v >> 8 );
        work[ 2 ] = ( byte ) ( v >> 16 );
        work[ 3 ] = ( byte ) ( v >> 24 );
        dis.write( work, 0, 4 );
        }

    /**
     * Write a long, 64-bits. like DataOutputStream.writeLong.
     *
     * @param v the long to write
     *
     * @throws IOException if write fails.
     */
    public final void writeLong( long v ) throws IOException
        {
        work[ 0 ] = ( byte ) v;
        work[ 1 ] = ( byte ) ( v >> 8 );
        work[ 2 ] = ( byte ) ( v >> 16 );
        work[ 3 ] = ( byte ) ( v >> 24 );
        work[ 4 ] = ( byte ) ( v >> 32 );
        work[ 5 ] = ( byte ) ( v >> 40 );
        work[ 6 ] = ( byte ) ( v >> 48 );
        work[ 7 ] = ( byte ) ( v >> 56 );
        dis.write( work, 0, 8 );
        }

    /**
     * Write short, 16-bits. Like DataOutputStream.writeShort. also acts as a writeUnsignedShort
     *
     * @param v the short you want written in little endian binary format
     *
     * @throws IOException if write fails.
     */
    public final void writeShort( int v ) throws IOException
        {
        work[ 0 ] = ( byte ) v;
        work[ 1 ] = ( byte ) ( v >> 8 );
        dis.write( work, 0, 2 );
        }

    /**
     * Write a string as a UTF counted string.
     *
     * @param s the string to write.
     *
     * @throws IOException if write fails.
     * @see java.io.DataOutput#writeUTF(java.lang.String)
     */
    public final void writeUTF( String s ) throws IOException
        {
        dis.writeUTF( s );
        }
    }// end LEDataOutputStream
