# !/usr/bin/env python

import ctypes
import os
import math
from numpy.ctypeslib import ndpointer, load_library
from numpy import empty, empty_like, require


def get_round(x):
    return int(x + 0.5) if x >= 0 else int(x - 0.5)


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
    xarray_type = ndpointer(float, flags='aligned, contiguous')
    oarray_type = ndpointer(float, ndim=2, flags='aligned, contiguous, writeable')
    carray_type = ctypes.c_char * len(compression)

    frameshift_samples = get_round(frameshift * float(fs) / 1000)
    numframes = int(math.ceil(numsamples / float(frameshift_samples)))

    x = require(x, float, ['CONTIGUOUS', 'ALIGNED'])
    result = empty((numframes, numchans), dtype=float)

    program_dirname = os.path.dirname(os.path.realpath(__file__))
    _libratemap = load_library('libratemap', program_dirname)
    _libratemap.ratemap.argtypes = (
        xarray_type, ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double,
        ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_char), oarray_type)

    _libratemap.ratemap.restype = None

    _libratemap.ratemap(x, ctypes.c_int(numsamples), ctypes.c_int(fs),
                        ctypes.c_double(lowcf), ctypes.c_double(highcf), ctypes.c_int(numchans),
                        ctypes.c_double(frameshift),
                        ctypes.c_double(ti), carray_type(*compression), result)

    return result


def gammatone(x, fs, cf, hrect=False):
    """
    Python wrapper for gammatone.c
    %  bm, env, instp, instf = gammatone_c(x, fs, cf, hrect)
    %
    %  x     - input signal
    %  fs    - sampling frequency (Hz)
    %  cf    - centre frequency of the filter (Hz)
    %  hrect - half-wave rectifying if hrect = True (default False)
    %
    %  bm    - basilar membrane displacement
    %  env   - instantaneous envelope
    %  instp - instantaneous phase (unwrapped radian)
    %  instf - instantaneous frequency (Hz)
    """

    numsamples = len(x)
    xarray_type = ndpointer(float, flags='aligned, contiguous')
    oarray_type = ndpointer(float, flags='aligned, contiguous, writeable')

    x = require(x, float, ['CONTIGUOUS', 'ALIGNED'])
    bm = empty_like(x)
    env = empty_like(x)
    instp = empty_like(x)
    instf = empty_like(x)

    program_dirname = os.path.dirname(os.path.realpath(__file__))
    _libratemap = load_library('libratemap', program_dirname)
    _libratemap.gammatone.argtypes = (
        xarray_type, ctypes.c_int,
        ctypes.c_int, ctypes.c_double, ctypes.c_int,
        oarray_type, oarray_type,
        oarray_type, oarray_type)

    _libratemap.gammatone.restype = None

    _libratemap.gammatone(x, ctypes.c_int(numsamples),
                          ctypes.c_int(fs), ctypes.c_double(cf), ctypes.c_int(hrect),
                          bm, env, instp, instf)

    return bm, env, instp, instf

