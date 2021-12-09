.PHONY : all
all : task1

task1 : task1.c
	gcc -Wall -g -o task1 task1.c -O3 -lm -lgsl -lgslcblas

clean : 
	touch *.c
