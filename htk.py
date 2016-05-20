"""htk
Reading and writing numpy arrays from/to HTK feature files
"""

import numpy as np
import struct

def htk_write(data_orig, filename, kind=9, byteorder='b', frame_period=0.01):
    """Write a numpy ndarray to an HTK format file"""
    data = data_orig.astype('f4')
    (nframes, ndims) = data.shape
    with open(filename,'w') as fid:
        fid.write(struct.pack('>iihh', nframes, frame_period*1e7, ndims*4, kind))
        # to avoid overflow when writing as 32 bit float
        #np.clip(data, np.finfo('float32').min, np.finfo('float32').max, data)
        np.clip(data, 0.0, 10.0, data)
        data[np.isnan(data)] = 0
        data[np.isinf(data)] = 0
        for row in data:
            fid.write(struct.pack('>%sf'%len(row), *row))

def htk_read(filename):
	"""Construct a numpy ndarray from an HTK format file"""
	with open(filename, 'r') as fid:
		(nframes, __, ndims, __) = struct.unpack('>iihh', fid.read(12))
		ndims /= 4
		data = np.empty((nframes, ndims))
		for i in xrange(len(data)):
			data[i, :] = struct.unpack('>%sf'%ndims, fid.read(ndims*4))
	return data


# Test code
if __name__ == '__main__':
	data = np.ones((32,12))*100;
	htk_write(data, 'test.rate32')
	x = htk_read('test.rate32')
	print(x)