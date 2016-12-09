CPPFLAGS = -I../S4 -I../S4/RNP -I../S4/kiss_fft

CLFLAGS = /c /GR /W3 /EHs /nologo /MD /D_USE_MATH_DEFINES /Dstrcasecmp=_stricmp /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /O2 /Oy- /DNDEBUG $(CPPFLAGS)

all:
	cl $(CLFLAGS) ..\S4\RNP\Eigensystems.cpp /FoC:Eigensystems.obj
	cl $(CLFLAGS) ..\S4\gsel.c /FoC:gsel.obj
	cl $(CLFLAGS) ..\S4\numalloc.c /FoC:numalloc.obj
	cl $(CLFLAGS) ..\S4\rcwa.cpp /FoC:rcwa.obj
	cl $(CLFLAGS) ..\S4\S4.cpp /FoC:S4.obj
	cl $(CLFLAGS) ..\S4\sort.c /FoC:sort.obj
	cl $(CLFLAGS) ..\S4\fmm\fft_iface.cpp /FoC:fft_iface.obj
	cl $(CLFLAGS) ..\S4\fmm\fmm_common.cpp /FoC:fmm_common.obj
	cl $(CLFLAGS) ..\S4\fmm\fmm_closed.cpp /FoC:fmm_closed.obj
	cl $(CLFLAGS) ..\S4\fmm\fmm_experimental.cpp /FoC:fmm_experimental.obj
	cl $(CLFLAGS) ..\S4\fmm\fmm_FFT.cpp /FoC:fmm_FFT.obj
	cl $(CLFLAGS) ..\S4\fmm\fmm_kottke.cpp /FoC:fmm_kottke.obj
	cl $(CLFLAGS) ..\S4\fmm\fmm_PolBasisJones.cpp /FoC:fmm_PolBasisJones.obj
	cl $(CLFLAGS) ..\S4\fmm\fmm_PolBasisNV.cpp /FoC:fmm_PolBasisNV.obj
	cl $(CLFLAGS) ..\S4\fmm\fmm_PolBasisVL.cpp /FoC:fmm_PolBasisVL.obj
	cl $(CLFLAGS) ..\S4\pattern\pattern.c /FoC:pattern.obj
	cl $(CLFLAGS) ..\S4\pattern\intersection.c /FoC:intersection.obj
	cl $(CLFLAGS) ..\S4\pattern\predicates.c /FoC:predicates.obj
	cl $(CLFLAGS) ..\S4\kiss_fft\kiss_fft.c /FoC:kiss_fft.obj
	cl $(CLFLAGS) ..\S4\kiss_fft\tools\kiss_fftnd.c /FoC:kiss_fftnd.obj
	lib /out:S4.lib *.obj

clean:
	del *.obj S4.lib
