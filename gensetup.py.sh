#!/bin/bash

OBJDIR="$1"
LIBFILE="$2"

cat <<SETUPPY > setup.py
from distutils.core import setup, Extension
import os

os.environ['CC'] = '/usr/bin/gcc'
os.environ['LDSHARED'] = '/usr/bin/gcc -shared'

S4module = Extension('S42',
	sources = [
		'S4/main_python.c'
	],
	libraries = [
		'S4',
		'stdc++',
                'openblas',
                'gfortran',
                'cholmod',
                'amd',
                'colamd',
                'camd',
                'ccolamd',
                'fftw3' 
	],
	library_dirs = ['$OBJDIR'],
	extra_link_args = [
		'$LIBFILE'
	], 
	extra_compile_args = ["-std=c99"],
)

setup(name = 'S42',
	version = '1.1',
	description = 'Stanford Stratified Structure Solver (S4): Fourier Modal Method',
	ext_modules = [S4module]
)
SETUPPY
