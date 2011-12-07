LIBS = /d/dev/libs/output/lapack-3.2.2-mingw32-vanilla/lapack_mingw32.a /d/dev/libs/output/lapack-3.2.2-mingw32-vanilla/blas_mingw32.a -lgfortran

PTHREAD_FLAGS =
PTHREAD_LIB = -lpthread -lwsock32

CFLAGS = -ggdb -O0

bad:
	g++ $(CFLAGS)        main.cpp $(PTHREAD_FLAGS) $(PTHREAD_LIB) $(LIBS)
good:
	g++ $(CFLAGS) -DGOOD main.cpp $(PTHREAD_FLAGS) $(PTHREAD_LIB) $(LIBS)
