.PHONY : all
all : task2

task2 : task2.c
	gcc -Wall -g -o task2 task2.c -lm -lgsl -lgslcblas

clean : 
	touch *.c
