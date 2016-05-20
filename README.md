# ratemap.py
A Python ctypes wrapper around Ning's ratemap implementation.

# Installation. 
Just compile c code with: cc -fPIC -shared -o libratemap.so ratemap.c

The python code uses the 'docopt' module which might need installing with 'pip install docopt'

# Usage:

Can be used to compute a single ratemap and save result to an HTK-format file, e.g.
./ratemap.py -o output.rate32 input.wav

Can be used to calculate and display a ratemap
./ratemap.py --display input.wav

It can process a 'script' file containing input/output filename pairs ala HTK
./ratemap.py -S scriptfile

Defaults to 32 channels form 50 to 3500 Hz. 

For help ,
./ratemap --help

# Note to drop it into the current HTK-python framework,:

./train.py -x ./ratemap.py -xp '--cuberoot -S' -T 0 -p data/spanish_phones_class.hed.orig -w data/spanish_wordlist.txt -a ../ES_wordset -d data/spanish_dictionary.txt -n "{1: 1000, 2:500}" -r 1
./adapt_cmllr.py  -x ./ratemap.py -xp '--cuberoot -S' -T 1 -n "{3: 1000, 4:500}" -r 1
./test.py  -x ./ratemap.py -xp '--cuberoot -S' -T 1 -n "{1: 1000, 2:500}" -r 1

# Issues:

- The gammatone filter can be unstable during low energy regions (?) Calculations overflow and produce nan's.  These are set to zero before writing to HTK. 
