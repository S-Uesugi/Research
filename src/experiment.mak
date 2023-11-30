CXX = g++
SRC = experiment.cpp kannan_tool.cpp weighted_LLL.cpp
OUT = -o experiment.out
CFLAGS = -pthread -march=native
LDFLAGS = -lntl -lgmp -lm 

all:
	${CXX} ${CFLAGS} ${OUT} ${SRC} ${LDFLAGS}