#!/bin/bash

OBJDIR="$1"
LIBFILE="$2"

cat <<SETUPPY > setup.py
from distutils.core import setup, Extension

S4module = Extension('S4',
	sources = [
		'S4/main_python.c'
	],
	libraries = [
		'S4',
		'stdc++'
	],
	library_dirs = ['$OBJDIR'],
	extra_link_args = [
		'$LIBFILE'
	]
)

setup(name = 'S4',
	version = '1.1',
	description = 'Stanford Stratified Structure Solver (S4): Fourier Modal Method',
	ext_modules = [S4module]
)
SETUPPY
