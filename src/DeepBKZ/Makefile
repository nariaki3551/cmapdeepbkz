CC	= g++
SRC	= main.cpp
# BOOST_DIR = ../boost_1_75_0
# BOOST_DIR = /home/tateiwa/usr/source/boost_1_80_0
BOOST_DIR = /home/tateiwa/usr/source/boost_1_75_0
CFLAGS	= -std=c++17 -march=native -Ofast -funroll-loops -flto -fomit-frame-pointer -fprefetch-loop-arrays -msse4 -fpermissive
LDFLAGS	= -lgmp -lmpfr -lgsl -lntl -fopenmp
INCL	= -I. -I${BOOST_DIR}

all:
	${CC} ${CFLAGS} ${SRC} ${INCL} ${LDFLAGS}


