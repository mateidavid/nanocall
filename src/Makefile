SHELL := /bin/bash

OPT_FLAG = -O0
CXXFLAGS = -Wall -Wextra -pedantic -g3 -ggdb -fno-eliminate-unused-debug-types ${OPT_FLAG}
M_CXXFLAGS = -std=c++11
CPPFLAGS = -I .
HDF_LIB = -lhdf5

TARGETS = compute-state-transitions compute-scaled-pore-model run-fwbw

.PHONY: all test clean

all: ${TARGETS}

clean:
	rm -rf ${TARGETS}

compute-state-transitions: compute-state-transitions.cpp logsum.cpp
	${CXX} ${M_CXXFLAGS} ${CXXFLAGS} ${CPPFLAGS} ${LDFLAGS} $^ -o $@

compute-scaled-pore-model: compute-scaled-pore-model.cpp
	${CXX} ${M_CXXFLAGS} ${CXXFLAGS} ${CPPFLAGS} -I ../fast5/src ${LDFLAGS} $^ -o $@ ${HDF_LIB}

run-fwbw: run-fwbw.cpp logsum.cpp
	${CXX} ${M_CXXFLAGS} ${CXXFLAGS} ${CPPFLAGS} ${LDFLAGS} $^ -o $@ -lz
