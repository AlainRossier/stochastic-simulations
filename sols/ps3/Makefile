INC	:= -I. -I..
LIB	:= -lm

all:	ps3

mcqmc06:	ps3.cpp Makefile \
		mlmc_test.cpp mlmc.cpp poissinv.h
		g++ -std=c++11 ps3.cpp -o ps3.o $(INC) $(LIB)

clean:
	rm -f ps3
