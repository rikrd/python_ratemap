from setuptools import setup, find_packages, Extension

libratemap = Extension('ratemap/libratemap',
                       sources=['ratemap/ratemap.c'])

setup(name='Ratemap',
      version='1.0',
      description='A Python ctypes wrapper around Ning''s ratemap implementation.',
      ext_modules=[libratemap],
      packages=find_packages(),
      install_requires=[i.strip() for i in open("requirements.txt").readlines()])
