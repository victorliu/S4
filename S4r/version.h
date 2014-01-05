#ifndef S4R_VERSION_H_INCLUDED
#define S4R_VERSION_H_INCLUDED

#define S4R_VERSION_MAJOR	"1"
#define S4R_VERSION_MINOR	"0"
#define S4R_VERSION_NUM		100
#define S4R_VERSION_RELEASE	"0"

#define S4R_VERSION	"S4r " S4R_VERSION_MAJOR "." S4R_VERSION_MINOR
#define S4R_RELEASE	S4R_VERSION "." S4R_VERSION_RELEASE
#define S4R_COPYRIGHT	S4R_RELEASE "  Copyright (C) 2013 Stanford University"
#define S4R_AUTHORS	"Victor Liu (vkl@alumni.stanford.edu)"

/* Reset the program name (default is "lua") */
#if defined(LUA_PROGNAME)
#undef LUA_PROGNAME
#define LUA_PROGNAME "S4r"
#endif

#endif /* S4R_VERSION_H_INCLUDED */
