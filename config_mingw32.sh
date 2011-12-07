#!/bin/sh

./configure --with-blas=/d/dev/libs/output/lapack-3.2.2-mingw32-vanilla/blas_mingw32.a --with-lapack=/d/dev/libs/output/lapack-3.2.2-mingw32-vanilla/lapack_mingw32.a CPPFLAGS="-I'/c/Progra~2/Lua/5.1/include' -I/d/dev/libs/output/fftw3 -I/d/dev/libs/output/numerical_mingw32/cholmod" LDFLAGS="-static -L/d/dev/libs/output/pthreads -L/d/dev/libs/output/lua-5.1.4/i386 -L/d/dev/libs/output/fftw3 -L/d/dev/libs/output/numerical_mingw32"
