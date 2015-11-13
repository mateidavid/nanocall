SHELL := /bin/bash

HDF_ROOT = /usr

OPT_FLAG = -O0
CXXFLAGS = -Wall -Wextra -pedantic -g3 -ggdb -fno-eliminate-unused-debug-types ${OPT_FLAG}
M_CXXFLAGS = -std=c++11
CPPFLAGS = -isystem ${HDF_ROOT}/include -I fast5/src -I tclap/include

TARGETS = compute-state-transitions compute-scaled-pore-model run-fwbw run-viterbi nanocall

.PHONY: all test clean

all: ${TARGETS}

clean:
	rm -rf ${TARGETS}

compute-state-transitions: compute-state-transitions.cpp logsum.cpp
	${CXX} ${M_CXXFLAGS} ${CXXFLAGS} ${CPPFLAGS} $^ -o $@ ${LDFLAGS}

compute-scaled-pore-model: compute-scaled-pore-model.cpp
	${CXX} ${M_CXXFLAGS} ${CXXFLAGS} ${CPPFLAGS} $^ -o $@ ${LDFLAGS} -L ${HDF_ROOT}/lib -lhdf5

run-fwbw: run-fwbw.cpp logsum.cpp
	${CXX} ${M_CXXFLAGS} ${CXXFLAGS} ${CPPFLAGS} $^ -o $@ ${LDFLAGS} -lz

run-viterbi: run-viterbi.cpp logsum.cpp
	${CXX} ${M_CXXFLAGS} ${CXXFLAGS} ${CPPFLAGS} $^ -o $@ ${LDFLAGS} -lz

nanocall: nanocall.cpp logsum.cpp Builtin_Modules.cpp
	${CXX} ${M_CXXFLAGS} ${CXXFLAGS} ${CPPFLAGS} $^ -o $@ ${LDFLAGS} -L ${HDF_ROOT}/lib -lhdf5 -lz
