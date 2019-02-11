#!/bin/bash

OBJDIR="$1"
LIBFILE="$2"
LIBS="$3"
x1="$4"
x2="$5"

cat <<SETUPPY > setup.py
from distutils.core import setup, Extension

libs = ['S4', 'stdc++']
libs.extend( [lib[2::] for lib in '$LIBS'.split()])

S4module = Extension('S4',
	sources = [
		'S4/main_python.c'
	],
	libraries = libs,
	library_dirs = ['$OBJDIR', $x2],
	extra_link_args = [
		'$LIBFILE', $x1
	],
	extra_compile_args=['-std=gnu99'] 
)

setup(name = 'S4',
	version = '1.1',
	description = 'Stanford Stratified Structure Solver (S4): Fourier Modal Method',
	ext_modules = [S4module]
)
SETUPPY
