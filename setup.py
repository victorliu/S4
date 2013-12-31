from distutils.core import setup, Extension
from distutils import ccompiler

S4module = Extension('S4',
	sources = [
		'src/main_python.c'
	],
	libraries = ['S4', 'stdc++'],
	library_dirs = ['./build'],
)

setup(name = 'S4',
	version = '1.0',
	description = 'Stanford Stratified Structure Solver (S4): Fourier Modal Method',
	ext_modules = [S4module]
)
