# Makefile
# Copyright (C) 2013 Michael Rose

PROJECT = ppmrose
CXX = g++
CXXFLAGS = -Wall -O3 -Isrc
TLIB = lib$(PROJECT).a
AR = ar
RM = rm -f
TARGET = ppmunwarp ppmwhitebalance


all: $(TARGET)

ppmunwarp: src/ppmunwarp.o $(TLIB)
	$(CXX) $(CXXFLAGS) $^ -o $@

ppmwhitebalance: src/ppmwhitebalance.o $(TLIB)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TLIB): src/ppmroselib.o
	$(AR) crs $@ $^

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGET) $(TLIB) src/*.o
