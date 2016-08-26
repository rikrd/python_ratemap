import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.io.wavfile as wavfile
import numpy as np
import ratemap

# Test code
if __name__ == '__main__':
    rate, x = wavfile.read('example.wav')
    y = ratemap.ratemap(x, rate, numchans=32)
    plt.imshow(np.flipud(y), cmap=m.Greys)
    plt.show()
