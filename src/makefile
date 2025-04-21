# CFLAGS for CC
CFLAGS=-g -O3 -std=c++17 -Wextra -fPIC# -stdlib=libc++
FCFLAGS=-w -c -fPIC -fno-automatic -std=legacy -ffree-line-length-none -O3

# Compiler
FF=gfortran $(FCFLAGS)
CCo=g++ -c $(CFLAGS)


# set the build directory
ifndef BUILDDIR 
	BUILDDIR=../build
endif

all: obj lib

obj:
	$(CCo) -o $(BUILDDIR)/t96.o t96.cc

lib: 
	$(CC) -fPIC $(BUILDDIR)/t96.o -shared -o ../lib/libt96.so
