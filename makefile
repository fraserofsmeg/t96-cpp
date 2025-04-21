
ifndef BUILDDIR 
	BUILDDIR=$(shell pwd)/build
endif

all: lib

lib: 
	mkdir -p $(BUILDDIR)
	mkdir -p lib
	cd src; make all

test:
	cd tests; make all


clean:
	-rm -rf build/
	-rm -rf lib/