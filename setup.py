from distutils.core import setup, Extension

ratemap = Extension('libratemap', sources=['ratemap.c'])

setup(name='Ratemap',
      version='1.0',
      description='A Python ctypes wrapper around Ning''s ratemap implementation.',
      ext_modules=[ratemap],
      py_modules=['ratemap', 'htk'],
      scripts=['ratemap.py'])
