# !/usr/bin/env python

import ctypes
import os
from numpy.ctypeslib import ndpointer


def ratemap(x, fs, lowcf=50.0, highcf=3500.0, numchans=32, frameshift=10.0,
            ti=8.0, compression='cuberoot'):
    """
    Python wrapper for ratemap.c
    %  ratemap = ratemap(x,fs,lowcf,highcf,numchans,frameshift,ti,compression)
    %
    %  x           input signal
    %  fs          sampling frequency in Hz
    %  lowcf       centre frequency of lowest filter in Hz (50)
    %  highcf      centre frequency of highest filter in Hz (3500)
    %  numchans    number of channels in filterbank (32)
    %  frameshift  interval between successive frames in ms (10)
    %  ti          temporal integration in ms (8)
    %  compression type of compression ['cuberoot','log','none'] ('cuberoot')
    """

    numsamples = len(x)
    xarray_type = ctypes.c_double * numsamples
    carray_type = ctypes.c_char * len(compression)

    frameshift_samples = int((frameshift * fs / 1000) + 0.5);
    numframes = (int)(1.0 + numsamples / frameshift_samples);

    program_dirname = os.path.dirname(os.path.realpath(__file__))
    _libratemap = ctypes.CDLL(program_dirname + '/libratemap.so')
    _libratemap.ratemap.argtypes = (
    ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double,
    ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_char))

    _libratemap.ratemap.restype = ndpointer(dtype=ctypes.c_double, shape=(numframes, numchans))

    result = _libratemap.ratemap(xarray_type(*x), ctypes.c_int(numsamples), ctypes.c_int(fs),
                                 ctypes.c_double(lowcf), ctypes.c_double(highcf), ctypes.c_int(numchans),
                                 ctypes.c_double(frameshift),
                                 ctypes.c_double(ti), carray_type(*compression))

    return result
