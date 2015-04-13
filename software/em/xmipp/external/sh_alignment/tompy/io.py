#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy

def read(filename):
    """Read EM file. Now only support read the type float32 on little-endian machines.

    @param filename: file name to read.

    @return The data from the EM file in ndarray
    """
    f = open(filename, 'r')
    try:
        dt_header = np.dtype('int32')
        header = np.fromfile(f, dt_header, 128)
        x = header[1]
        y = header[2]
        z = header[3]

        # read the data type
        dt = int(hex(header[0])[2])
        default_type = False

        if dt == 1: # byte
            raise Exception("Data type not supported yet!")
        elif dt == 2: # short
            dt_data = np.dtype('<i2')
        elif dt == 4: # long
            dt_data = np.dtype('<i4')
        elif dt == 5: # float32
            dt_data = np.dtype('<f4') # little-endian, float32
            default_type = True
        elif dt == 8: # float complex
            raise Exception("Data type not supported yet!")
        elif dt == 9: # double
            dt_data = np.dtype('<f8')
        elif dt == 10: # double complex
            raise Exception("Data type not supported yet!")
        else:
            raise Exception("Data type not supported yet!")

        v = np.fromfile(f, dt_data, x*y*z)
    finally:
        f.close()

    if default_type:
        volume = v.reshape((x, y, z), order='F') # fortran-order array
    else: # if the input data is not the default type, convert
        volume = np.array(v.reshape((x, y, z), order='F'), dtype='float32') # fortran-order array
    return volume


def write(filename, data):
    """Write EM file. Now only support written in type float32 on little-endian machines.

    @param filename: file name.
    @param data: data to write.
    """
    if data.dtype != np.dtype('float32'): # if the origin data type is not float32, convert
        data = np.array(data, dtype='float32')
    
    header = np.zeros(128, dtype='int32')
    header[0] = 83886086 # '0x5000006', TODO: hard-coded, to be changed!

    if len(data.shape) == 3:
        header[1:4] = data.shape
    elif len(data.shape) == 2:
        header[1:3] = data.shape
        header[3] = 1
    else:
        raise Exception("Input data shape invalid!")

    f = open(filename, 'wb')
    try:
        f.write(header.tostring())
        f.write(data.tostring(order='F')) # fortran-order array
    finally:
        f.close()


def n2v(data):
    """Transfer Numpy array into a Pytom volume.
    Note the data will be shared between Numpy and Pytom!

    @param data: data to convert.
    """
    try:
        from pytom_volume import vol
        from pytom_numpy import npy2vol
    except:
        raise ImportError("Pytom library is not installed or set properly!")

    if data.dtype != np.dtype("float32"):
        data = np.array(data, dtype="float32")
    
    if len(data.shape) == 3:
        if np.isfortran(data):
            # Fortran order
            v = npy2vol(data, 3)
        else:
            vv = np.asfortranarray(data)
            v = npy2vol(vv, 3)
    elif len(data.shape) == 2:
        if np.isfortran(data):
            v = npy2vol(data, 2)
        else:
            vv = np.asfortranarray(data)
            v = npy2vol(vv, 2)
    else:
        raise Exception("Data shape invalid!")

    return v

