#!/usr/bin/env python
""" ratemap.py 
Compute a ratemap from a wavefile. Output can be either displayed or saved to HTK file.

Usage:
  ratemap.py [-l <hz>][-h <hz>] [-n <nchans>] [-f <ms>] [-t <ms>] [--cuberoot|--log] [--display] [-o <file>] <wavfile_in>
  ratemap.py [-l <hz>][-h <hz>] [-n <nchans>] [-f <ms>] [-t <ms>] [--cuberoot|--log] [--display] -S <script_file>
  ratemap.py --help

Option:
  <wavefile_in>                        the wave file to process
  -l <hz>, --low=<hz>                  centre frequency of lowest filter in Hz [default: 50]
  -h <hz>, --high=<hz>                 centre frequency of highest filter in Hz [default: 3500]
  -n <chans>, --nchans=<chans>         number of ratemap channels [default: 32]
  -f <ms>, --frameshift=<ms>           interval between successive frames in ms [default: 10]
  -t <ms>, --temp_int=<ms>             temporal integration in ms [default: 8]
  -o <file>, --outfile==<file>         write the ratemap to an HTK format file
  -S <scriptfile>                      process a list of input/output file pairs
  --cuberoot        apply cuberoot compression 
  --log             apply log compression 
  --display         display the ratemap
  --help            print this help screen
"""

import ctypes
import os
import shlex
from functools import partial

import docopt
import numpy as np
import scipy.io.wavfile as wavfile
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

    frameshift_samples = int((frameshift*fs/1000)+0.5);
    numframes = (int)(1.0 + numsamples / frameshift_samples);
 
    program_dirname = os.path.dirname(os.path.realpath(__file__))
    _libratemap = ctypes.CDLL(program_dirname + '/libratemap.so')
    _libratemap.ratemap.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double, 
    ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_char))

    _libratemap.ratemap.restype = ndpointer(dtype=ctypes.c_double, shape=(numframes, numchans))

    result = _libratemap.ratemap(xarray_type(*x), ctypes.c_int(numsamples), ctypes.c_int(fs),
        ctypes.c_double(lowcf), ctypes.c_double(highcf), ctypes.c_int(numchans), ctypes.c_double(frameshift), 
        ctypes.c_double(ti), carray_type(*compression))

    
    return result


def process_one_file(make_features, infile, outfile, display):
    """Process a single file. Called in a loop if -S is used"""
    rate, x = wavfile.read(infile)
    y = make_features(x, rate)
    if display:
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        plt.imshow(np.flipud(np.transpose(y)), cmap = cm.Greys)
        plt.show()
    if outfile:
        import htk
        htk.htk_write(y, outfile)


if __name__ == '__main__':

    arguments = docopt.docopt(__doc__)

    compression='none'
    if arguments['--log']:
        compression = 'log'
    elif arguments['--cuberoot']:
        compression = 'cuberoot'

    lowcf = float(arguments['--low'])
    highcf = float(arguments['--high'])
    numchans = int(arguments['--nchans'])
    frameshift = float(arguments['--frameshift'])
    ti = float(arguments['--temp_int'])

    my_ratemap = partial(ratemap, lowcf=lowcf, highcf=highcf, numchans=numchans, 
        frameshift=frameshift, ti=ti, compression=compression)

    display = arguments['--display']
    outfile = arguments['--outfile']    

    if not arguments['-S']:
        infile = arguments['<wavfile_in>']
        process_one_file(my_ratemap, infile, outfile, display)
    else:
        scriptfile = arguments['-S']
        x = [tuple(shlex.split(line)) for line in open(scriptfile)]
        for (infile, outfile) in x:
            process_one_file(my_ratemap, infile, outfile, display)
