## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

## Paths to include
#AM_CPPFLAGS = -Wall -fpic -O2 -std=c++11 -pthread
AM_CPPFLAGS = -Wall -fpic -g -O0 -std=c++11 -pthread

noinst_HEADERS = arguments.h  bed.h  bg.h  Constant.h  Error.h  find_tss.h  find_tss_hk.h params.h  smooth.h  tabfile.h pException.h

bin_PROGRAMS = find_tss find_tss_hk

find_tss_SOURCES = find_tss.cpp find_tss_func.cpp arguments.cpp bg.cpp tabfile.cpp smooth.cpp bed.cpp

find_tss_hk_SOURCES = find_tss_hk.cpp find_tss_func_hk.cpp arguments.cpp bg.cpp tabfile.cpp smooth.cpp bed.cpp params.cpp Error.cpp
