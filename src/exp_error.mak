CXX = g++
SRC = exp_error.cpp kannan_tool.cpp weighted_LLL.cpp
OUT = -o exp_error.out
CFLAGS = -pthread -march=native
LDFLAGS = -lntl -lgmp -lm 

all:
	${CXX} ${CFLAGS} ${OUT} ${SRC} ${LDFLAGS}