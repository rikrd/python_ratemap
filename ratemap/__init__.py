# !/usr/bin/env python

import ctypes
import os
import math
from numpy.ctypeslib import ndpointer, load_library
from numpy import empty, empty_like, require
import numpy as np

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


def hz_to_erb_rate(x):
    return 21.4 * np.log10(4.37e-3 * x + 1.)


def erb_rate_to_hz(x):
    return (10**(x/21.4) - 1.) / 4.37e-3


def make_erb_cfs(lowcf, highcf, numchans):
    return erb_rate_to_hz(np.linspace(hz_to_erb_rate(lowcf), hz_to_erb_rate(highcf), numchans))


def hamming_win(n):
    return 0.54 - 0.46 * np.cos(2 * np.pi / (n - 1.) * np.arange(n, dtype=float))


def frames_to_samples(n, fr, fs):
    """ Convert from frame n at frame rate fr ms, to samples at samplrate Hz

    :param n: frame index
    :param fr: frame rate in ms
    :param fs: samplerate in Hz
    :return: sample index
    """

    return round((fs * fr * (n - 1)) / 1000) + 1


def synthesize_ratemap(x, fs, ratemap,
                       lowcf=50.0, highcf=3500.0,
                       numchans=32, frameshift=10.0, winlength=25.0,
                       compression='cuberoot'):
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

    cf = make_erb_cfs(lowcf, highcf, numchans)

    # TODO: Check that this framer is consistent to the one in gammatone.c
    winsize = int(round(fs * winlength / 1000.))  # in samples
    win = hamming_win(winsize)
    win4 = round(winsize / 4.)
    frameshift_samples = get_round(fs * frameshift / 1000.)
    num_frames = int(math.floor(numsamples / frameshift_samples) - math.ceil(winsize / frameshift_samples) + 1)

    s = np.zeros_like(x)

    if compression == 'log':
        ratemap_decomp = 10. ** ratemap

    elif compression == 'cuberoot':
        ratemap_decomp = ratemap ** 3

    else:
        raise ValueError("Only compression types 'log' and 'cuberoot' are supported.")

    for c in range(numchans):
        bm, _, _, _ = gammatone(x, fs, cf[c])

        # time reverse and do it again
        bm = bm[::-1]
        bm, _, _, _ = gammatone(bm, fs, cf[c])
        bm = bm[::-1]

        # overlap add in frames
        for frame in range(num_frames):
            idx1 = int(max(0, frames_to_samples(frame, frameshift, fs) - win4 - 1))
            idx2 = int(min(numsamples-1, idx1 + winsize))
            framelen = idx2 - idx1

            # If we use the energy of the window size here, we obtain worse results
            energy_src = np.sqrt(sum((bm[idx1:idx2]) ** 2))
            energy_tgt = ratemap_decomp[frame, c]

            s[idx1:idx2] = s[idx1:idx2] + bm[idx1:idx2] * win[:framelen] * energy_tgt / energy_src

    return s


def ratemap_for_synthesis(x, fs,
                          lowcf=50.0, highcf=3500.0,
                          numchans=32, frameshift=10.0, winlength=25.0,
                          compression='cuberoot'):
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

    cf = make_erb_cfs(lowcf, highcf, numchans)

    # TODO: Check that this framer is consistent to the one in gammatone.c
    winsize = int(round(fs * winlength / 1000.))  # in samples

    # TODO: Check why use ones for analysis and Hamming for synthesis
    win = hamming_win(winsize)
    win = np.ones(winsize)
    win4 = round(winsize / 4.)
    frameshift_samples = get_round(fs * frameshift / 1000.)
    num_frames = int(math.floor(numsamples / frameshift_samples) - math.ceil(winsize / frameshift_samples) + 1)

    ratemap = np.zeros((num_frames, numchans))

    for c in range(numchans):
        bm, _, _, _ = gammatone(x, fs, cf[c])

        # time reverse and do it again
        bm = bm[::-1]
        bm, _, _, _ = gammatone(bm, fs, cf[c])
        bm = bm[::-1]

        # overlap add in frames
        for frame in range(num_frames):
            idx1 = int(max(0, frames_to_samples(frame, frameshift, fs) - win4 - 1))
            idx2 = int(min(numsamples-1, idx1 + winsize))
            framelen = idx2 - idx1

            energy = np.sqrt(sum((bm[idx1:idx2] * win[:framelen]) ** 2))
            ratemap[frame, c] = energy

    if compression == 'log':
        ratemap_comp = np.log10(ratemap)

    elif compression == 'cuberoot':
        ratemap_comp = np.power(ratemap, 1/3.)

    else:
        raise ValueError("Only compression types 'log' and 'cuberoot' are supported.")

    return ratemap_comp
