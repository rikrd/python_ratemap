from setuptools import setup, find_packages, Extension

libratemap = Extension('ratemap/libratemap',
                       sources=['ratemap/ratemap.c',
                                'ratemap/gammatone.c'])

setup(name='Ratemap',
      version='1.1',
      author='Ricard Marxer',
      description='A Python ctypes wrapper around Ning''s ratemap implementation.',
      ext_modules=[libratemap],
      packages=find_packages(),
      install_requires=[i.strip() for i in open("requirements.txt").readlines()])
