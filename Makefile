# ============================================================================ 
# Name        : Makefile
# Author      : Ricardo Serrano
# Version     :
# Copyright   : This is CopyRighted to Ricardo Serrano
# Description : Makefile for Explicit FEM with MPI in Fortran
# ============================================================================

.PHONY: all clean

.SUFFIXES: .f90

OS := $(shell uname)

GCCVERSION = $(shell gcc --version | grep ^gcc | sed 's/^.* //g' | cut -c 1-3)

FC=mpif90 -fno-range-check -O3 -Wall -cpp -Wno-unused-variable 
CPPC=g++ -O3 -Wall

ifeq ($(OS), Darwin)
LIBS=-L./libs -framework accelerate -lpthread -lparmetis -lmetis -lstdc++
endif
ifeq ($(OS), Linux)
LIBS=-L./libs -lacml -lacml_mv -lpthread -lparmetis -lmetis -lstdc++ -lc
endif

FLAGS=
FLAGS+=$(INCS) $(LIBS)

SOURCES=utilities.f90 io_utils.f90 formulation.f90 pwaves.f90 model.f90 mesh_utils.f90 communicator.f90 explicit.f90 quadtree.f90 visualization.f90 test.f90 main.f90 
OBJECTS=$(SOURCES:.f90=.o) io.o

.f90.o:
	$(FC) $(INCS) -c $< -o $@

all debug: exfem

exfem: util_temps.f90 comm_temps.f90 $(OBJECTS)
	$(FC) $(INCS) -o $@ $(OBJECTS) $(LIBS)

comm_temps.f90: comm_temps.m4 
	m4 comm_temps.m4 > comm_temps.f90
util_temps.f90: util_temps.m4
	m4 util_temps.m4 > util_temps.f90
io.o: io.cpp
	$(CPPC) -c $< 
clean:
	rm -f exfem *.mod *.o comm_temps.f90 util_temps.f90 exfem_run.*

debug:	FC=mpif90 -fno-range-check -g3 -Wall -Wno-unused-variable -cpp

debug:	CPPC=g++ -g3 -Wall

