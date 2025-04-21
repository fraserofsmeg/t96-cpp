# CFLAGS for CC
CFLAGS=--g -O3 -Wextra -fPIC -I../src
CXXFLAGS=-g -O3 -std=c++17 -Wextra -fPIC -I../src
FCFLAGS=-w -c -fPIC -fno-automatic -std=legacy -ffree-line-length-none -O3


OS=$(shell uname -s)
CC=gcc $(CFLAGS)
CXX=g++ $(CXXFLAGS)
FF=gfortran $(FCFLAGS)
LDFLAGS=-lm -Wl,-rpath='$$ORIGIN/../lib' -L ../lib -lt89 


all:
	$(FF) -o t89c.o T89c.f 
	$(CXX) -o testfortran testfortran.cc t89c.o $(LDFLAGS)
	LD_LIBRARY_PATH=../lib ./testfortran
	-rm -v t89c.o
	-rm -v testfortran
