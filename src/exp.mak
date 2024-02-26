CXX = g++
SRC = exp.cpp kannan_tool.cpp weighted_LLL.cpp
OUT = -o exp.out
CFLAGS = -pthread -march=native
LDFLAGS = -lntl -lgmp -lm 

all:
	${CXX} ${CFLAGS} ${OUT} ${SRC} ${LDFLAGS}