# Makefile
# Copyright (C) 2013 Michael Rose

all: ppmunwarp ppmwhitebalance

ppmroselib.o:  ppmroselib.cc ppmroselib.h
	g++ -Wall -O3 -c ppmroselib.cc

ppmunwarp:  ppmunwarp.cc ppmroselib.o ppmroselib.h
	g++ -Wall -O3 -o ppmunwarp ppmroselib.o ppmunwarp.cc

ppmwhitebalance:  ppmwhitebalance.cc ppmroselib.o ppmroselib.h
	g++ -Wall -O3 -o ppmwhitebalance ppmroselib.o ppmwhitebalance.cc

clean:
	rm -f ppmroselib.o ppmunwarp ppmwhitebalance
