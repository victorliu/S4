ifeq ($(OS),Windows_NT)
    include Makefile.Msys2
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        include Makefile.Linux
    endif
    ifeq ($(UNAME_S),Darwin)
        include Makefile.Darwin
    endif
endif

include Makefile.common
