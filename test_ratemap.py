import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.io.wavfile as wavfile
import numpy as np
import ratemap


def test_ratemap():
    rate, x = wavfile.read('example.wav')
    y = ratemap.ratemap(x, rate, numchans=32)
    plt.imshow(y.T, cmap=cm.Greys, aspect='auto', interpolation='nearest', origin='lower')
    plt.show()


def test_synthesize_ratemap():
    rate, x = wavfile.read('example.wav')
    for compression in ['cuberoot', 'log']:
        y = ratemap.ratemap_for_synthesis(x, rate, numchans=32, compression=compression)
        x_hat = ratemap.synthesize_ratemap(x, rate, y, numchans=32, compression=compression)
        wavfile.write('example_resynth_{}.wav'.format(compression), rate, x_hat)


def test_gammatone():
    rate, x = wavfile.read('example.wav')
    bm, env, instp, instf = ratemap.gammatone(x, rate, 150.0)
    plt.plot(bm[:1000])
    plt.show()


if __name__ == '__main__':
    # test_gammatone()
    test_ratemap()
    test_synthesize_ratemap()
