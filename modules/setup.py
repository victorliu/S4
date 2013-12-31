from distutils.core import setup, Extension

module1 = Extension('FunctionSampler1D',
                    sources = ['py_function_sampler_1d.c','function_sampler_1d.c'])
module2 = Extension('FunctionSampler2D',
                    sources = ['py_function_sampler_2d.c','function_sampler_2d.c','predicates.c'])

setup (name = 'FunctionSampler1D',
       version = '1.0',
       description = 'Adaptive resolution sampler for (smooth/Lipschitz) functions in 1D',
       ext_modules = [module1])
setup (name = 'FunctionSampler2D',
       version = '1.0',
       description = 'Adaptive resolution sampler for (smooth/Lipschitz) functions in 2D',
       ext_modules = [module2])
