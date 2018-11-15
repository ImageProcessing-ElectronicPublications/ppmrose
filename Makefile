# Makefile
# Copyright (C) 2013 Michael Rose

all: ppmunwarp

ppmunwarp:  ppmunwarp.cc
	g++ -Wall -O3 -o ppmunwarp ppmunwarp.cc

clean:
	rm -f ppmunwarp
